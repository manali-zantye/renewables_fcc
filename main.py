import numpy as np
import pandas as pd
import os
import sqlite3

#read power plant data
plant_data = pd.read_csv('data_files/power_plant_data.csv', index_col=0)

#select a plant for demo
sel_plant = plant_data.loc[148]

##%%%Algorithm Step 0: Get time-series data for electricity price, solar and wind capacity factors%%###
all_weather = pd.read_csv('data_files/weather_station_metadata.csv', index_col=0) #data for all weather stations

def haversine_distance(lat1, lon1, lat2, lon2):
   r = 6371
   phi1 = np.radians(lat1)
   phi2 = np.radians(lat2)
   delta_phi = np.radians(lat2 - lat1)
   delta_lambda = np.radians(lon2 - lon1)
   a = np.sin(delta_phi / 2)**2 + np.cos(phi1) * np.cos(phi2) *   np.sin(delta_lambda / 2)**2
   res = r * (2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a)))
   return np.round(res, 2)

temp_weather = all_weather.copy() #get data for closest weather station
distances_km = []    
for row in temp_weather.iterrows():
    distances_km.append(
    haversine_distance(sel_plant.loc['Latitude'], sel_plant.loc['Longitude'], row[1].Latitude, row[1].Longitude)
    )    
temp_weather['dis_frompp'] = distances_km
closest_stat = temp_weather['dis_frompp'].idxmin()
xls = pd.ExcelFile('data_files/weather_data.xlsx')
df_weather = pd.read_excel(xls, str(int(closest_stat)), index_col=0)

df_price = pd.read_csv("data_files/electricity_price.csv", index_col=0)
agg_price = df_price.groupby(['Hour']).mean()

##%%%Algorithm Step 1: Calculate parameter vector, formulate scenario tree%%%##

vci = 1.5
vr = 12
vco = 25
df_weather.set_index(agg_price.index, inplace=True)

# #calc capacity factor of solar PV for each data point
df_weather['PV_capfac'] = (df_weather['GHI (W/m^2)']/1000)*0.9375

# #calc capacity factor of wind turbines for each data point

df_weather.loc[(df_weather['Wspd (m/s)'] >= vci) & (df_weather['Wspd (m/s)'] <= vr), 'W_capfac'] = (df_weather['Wspd (m/s)']**3 - vci**3)/(vr**3 - vci**3)
df_weather.loc[(df_weather['Wspd (m/s)'] <= vci) | (df_weather['Wspd (m/s)'] >= vco), 'W_capfac'] = 0
df_weather.loc[(df_weather['Wspd (m/s)'] >= vr) & (df_weather['Wspd (m/s)'] <= vco), 'W_capfac'] = 1
weather_write = df_weather[['PV_capfac','W_capfac']]
combined_write = pd.concat([agg_price, weather_write], axis=1)
combined_write.columns = ['price', 'PV_capfac','W_capfac']
combined_write.to_csv('parameters.csv')

#Write power plant nameplate capacity and emission intensity 
f = open('namepcap.gms', "w+")
f.write("%f" % (sel_plant.loc['Nameplate capacity (MW)']))
f.close()

f = open('co2eminten.gms', "w+")
f.write("%f" % (sel_plant.loc['Emission intensity (ton/MWh)']))
f.close()

data_scenarios = combined_write[['price', 'W_capfac', 'PV_capfac']]

#Formulate scenario tree

num_scenarios = len(data_scenarios.index)
num_processes = len(data_scenarios.columns)
prob = (1/num_scenarios)*np.ones((num_scenarios, num_processes))

#nodes demands
f= open("scen_red/nodes_demands.gms", "w+")
for i in range(num_processes):
        for j in range(num_scenarios):
                f.write("%d %s\r\n" % (i*num_scenarios + j + 2, str(round(data_scenarios.iloc[j,i],7)))) 
f.close()

#nodes probabilities
f = open("scen_red/nodes_probabilities.gms", "w+")
for i in range(num_processes):
        for j in range(num_scenarios):
                f.write("%d %s\r\n" % (i*num_scenarios + j + 2, str(prob[j,i]))) 
f.close()


##%%%Algorithm Step 2: Specify scenario reduction accuracy%%###
red_accuracy = 0.94
f = open("scen_red/SR_category.gms", "w+")
f.write("ScenRedParms('red_percentage') = %1.3f;" % (1 - red_accuracy))
f.close()

##%%%Algorithm Step 3: Perform scenario reduction and obtained reduced data%%###

os.system('cd scen_red;gams my_SR_code.gms u1 = 26281 u2 = 8760 u3 = 4;cd ..')

#functions to read reduced data 
def read_reduced_scenario_data(fpath, num_processes):
    results_data = np.genfromtxt(fpath, delimiter=',')
    num_reduced_scen = int(np.size(results_data)/(3*num_processes))

    r_nod = results_data[0:num_processes*num_reduced_scen].reshape(num_processes,num_reduced_scen)
    r_nod = np.transpose(r_nod)

    r_prob = results_data[num_processes*num_reduced_scen:2*num_processes*num_reduced_scen].reshape(num_processes,num_reduced_scen)
    r_prob = np.transpose(r_prob)

    r_dem = results_data[2*num_processes*num_reduced_scen:3*num_processes*num_reduced_scen].reshape(num_processes,num_reduced_scen)
    r_dem = np.transpose(r_dem)
    
    return r_nod,r_prob,r_dem,num_reduced_scen


def write_gams_file_main(num_scen, data_scenarios, time_periods):  
   
    f = open("red_opti/price.txt","w+")
    for j in range(num_scen):
        f.write("%d %s\n" % ((j+1),str(data_scenarios[j,0])))
    f.close()

    f = open("red_opti/cap_fac.txt","w+")
    for j in range(num_scen):
        f.write("w.%d %s\n" % ((j+1),str(data_scenarios[j,1])))
    f.close()
    

    f = open("red_opti/cap_fac.txt","a")
    for j in range(num_scen):
        f.write("sp.%d %s\n" % ((j+1),str(data_scenarios[j,2])))
    f.close()
    
    f = open("red_opti/solar_rad.txt","w+")
    for j in range(num_scen):
        f.write("%d %d\n" % ((j+1),0))
    f.close()

    f = open("red_opti/time_periods.txt","w+")
    for j in range(num_scen):
        f.write("%d %s\n" % ((j+1),str(time_periods[j])))
    f.close()
    
    return

#read reduced scenario data

num_time_periods = 24*365
num_processes = 3
scen_file = 'scen_red/results.csv' #read reduced scenario data
r_nod,r_prob,r_dem,num_reduced_scen = read_reduced_scenario_data(scen_file, num_processes)
r_time_periods = r_prob[:,0]*num_time_periods
write_gams_file_main(num_reduced_scen, r_dem, r_time_periods)

##%%%Algorithm Step 4: Use reduced scenario data to obtain optimal design%%%##

# Considering future cost case: write cost parameters
price_ren = 0.3
tax_carbon = 80
eor_price = 35

f = open("red_opti/price_pv.gms","w+")
f.write(str(price_ren) + "e6")
f.close()

f = open("red_opti/price_wt.gms","w+")
f.write(str(price_ren) + "e6")
f.close()

f = open("red_opti/tax_carbon.gms","w+")
f.write("%d" % (tax_carbon))
f.close()

f = open("red_opti/eor_price.gms","w+")
f.write("%d" % (eor_price))
f.close()

#copy power plant nameplate capacity and emission  intensity
os.system('scp namepcap.gms red_opti')
os.system('scp co2eminten.gms red_opti')

#run GAMS file for optimization
os.system('cd red_opti;gams main.gms u1 = ' + str(num_reduced_scen+1) + ';cd ..')

#read optimal design data
results_file = 'results.db'
database = 'red_opti/' + results_file
dat = sqlite3.connect(database)
scal_var_gams = pd.read_sql_query("SELECT * FROM scalarvariables", dat).set_index('name').level

results_df = pd.DataFrame({'y_capt': [scal_var_gams.loc['y_capt'].astype(int)], 
                           'sz_pv': [pd.read_sql_query("SELECT * FROM sz", dat).set_index('m').level.loc['sp'].astype(float)],
                           'sz_w': [pd.read_sql_query("SELECT * FROM sz", dat).set_index('m').level.loc['w'].astype(float)],
                           'npv': [scal_var_gams.loc['NPV'].astype(float)],
                           'sz_b': [scal_var_gams.loc['nc_b'].astype(float)]})


#write designs to fix_design.gms file
f = open("red_opti/fix_design.gms", "w+")
f.write("y_capt.fx = %d;\n" %(results_df['y_capt']))
f.write("sz.fx('w') = %f;\n" %(results_df['sz_w']))
f.write("sz.fx('sp') = %f;\n" %(results_df['sz_pv']))
f.write("nc_b.fx = %f;\n" %(results_df['sz_b']))
f.close()

##%%Algorithm Step 5: Obtain optimal operational decisions using original dataset%%##

#copy relevant files to the new folder
os.system('scp red_opti/fix_design.gms orig_opti')
os.system('scp namepcap.gms orig_opti')
os.system('scp co2eminten.gms orig_opti')
os.system('scp red_opti/price_pv.gms orig_opti')
os.system('scp red_opti/price_wt.gms orig_opti')
os.system('scp red_opti/tax_carbon.gms orig_opti')
os.system('scp red_opti/eor_price.gms orig_opti')

#load original, unreduced data

orig_data = data_scenarios.reset_index()
orig_data['Hour'] = orig_data['Hour'] + 1
orig_data = orig_data.set_index('Hour', drop=True)
orig_data['time_periods'] = np.ones(len(orig_data.index))
orig_data['solar_rad'] = np.zeros(len(orig_data.index))
orig_data['price'].to_csv('orig_opti/price.txt', header=None, sep=' ')
orig_data['time_periods'].to_csv('orig_opti/time_periods.txt', header=None, sep=' ')
orig_data['solar_rad'].to_csv('orig_opti/solar_rad.txt', header=None, sep=' ')


f = open('orig_opti/cap_fac.txt', 'w+')
for i in (orig_data.index):
        f.write('w.' + str(i) + ' ' + str(orig_data.loc[i,'W_capfac']) + '\n')
f.close()


f = open('orig_opti/cap_fac.txt', 'a')
for i in (orig_data.index):
        f.write('sp.' + str(i) + ' ' + str(orig_data.loc[i,'PV_capfac']) + '\n')
f.close()

#run GAMS file to obtain optimal operation
os.system('cd orig_opti;gams main.gms u1 = ' + str(8761) + ';cd ..')
