#
# This is the classic 1D Carbon box mode.
# for predicting atmospheric carbon dioxide
# concentrations
#
# This is based on the following paper:
#
# Rodhe, H. and Bjorkstrom, A. Some consequences of
# non-proportionality between fluxes and reservoir contents in
# natural systems, Tellus (1979), 31, 269-278
# https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.2153-3490.1979.tb00905.x
#
# Will Grey
# 2015-03-23
#

import matplotlib.pyplot as plt

PPM2GT = 2.123

def emissionsScenario(scenario):

    # Data from IPCC

  # These data are projected emissions for the years
  # 1990, 2000, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100
    
    if scenario=="A1":
        e=[7.10,7.97,10.88,12.64,14.48,15.35,16.38,16.00,15.73,15.18,14.30,13.49]
    if scenario=="A1G":
        e=[7.10,7.97,9.73,12.73,16.19,19.97,23.90,25.69,27.28,28.68,28.42,28.24]
    if scenario=="A1T":
        e= [7.10,7.97,9.38,10.26,12.38,12.65,12.26,11.38,9.87,8.02,6.26,4.32]
    if scenario=="A2":
        e= [7.10,7.97,9.58,12.25,14.72,16.07,17.43,19.16,20.89,23.22,26.15,29.09]
    if scenario=="B1":
        e= [7.10,7.97,9.28,10.63,11.11,11.72,11.29,9.74,8.18,6.70,5.32,4.23]
    if scenario=="B2":
        e= [7.10,7.97,8.78,9.05,9.90,10.69,11.01,11.49,11.62,12.15,12.79,13.32]

    fullEmissionScenario=[]

    # Interpolate between each of the ten years

    for i in range(11):
        for j in range(10):
             
             fullEmissionScenario.append(e[i]*(1-(j/10))+e[i+1]*(j/10))
        
    return  fullEmissionScenario



def calc_total_anthropogenic_emissions(e="Total"):

    # data from http://www.globalcarbonproject.org/carbonbudget
    # 1959-2017
    
    fossil_fuel_emissions=[2.45, 2.57, 2.58, 2.69, 2.83, 2.99, 3.13, 3.29,
                                 3.39, 3.57, 3.78, 4.05, 4.21 ,4.37, 4.61,
                                 4.62, 4.59, 4.86, 5.01, 5.07, 5.35, 5.29,
                                 5.13, 5.08, 5.06, 5.24, 5.40, 5.57, 5.71,
                                 5.92, 6.05, 6.05, 6.12, 6.06, 6.05, 6.15,
                                 6.28, 6.42, 6.53, 6.55, 6.53, 6.70, 6.86,
                                 6.96, 7.33, 7.69, 7.98, 8.27, 8.43, 8.70,
                                 8.60, 9.02, 9.38, 9.53, 9.61, 9.69, 9.68, 9.74, 9.87]

    land_use_change_emissions=[1.81, 1.67, 1.61, 1.57, 1.51, 1.46,
                                             1.42, 1.38, 1.35, 1.35, 1.36, 1.32, 1.28, 1.25,
                                             1.22, 1.18, 1.15, 1.13, 1.11, 1.07, 1.04, 1.09,
                                             1.10, 1.11, 1.15, 1.19, 1.20, 1.24, 1.28, 1.30,
                                             1.32, 1.32, 1.33, 1.35, 1.35, 1.34, 1.33, 1.31,
                                             1.78, 1.23, 1.20, 1.32, 1.18, 1.34, 1.39, 1.34,
                                             1.21, 1.25, 1.07, 1.13, 1.57, 1.42, 1.36, 1.60,
                                             1.54, 1.60, 1.62, 1.30, 1.39]

    total_anthropogenic_emissions=[]
    for i in range(len(fossil_fuel_emissions)):
        total_anthropogenic_emissions.append(fossil_fuel_emissions[i]+land_use_change_emissions[i])

    if e=="Total":
        return total_anthropogenic_emissions
    elif e=="FF":
        return fossil_fuel_emissions
    elif e=="LUC":
        return land_use_change_emissions
    
        

def box_model(initialAtmosCarbon, years, total_anthropogenic_emissions):        		

    # Time Step (dt)          
    dt=1

    # Initial values (Gt) of Carbon reservoir stores
    atmosCarbonStoreStart =initialAtmosCarbon*PPM2GT   
    surfaceOceanCarbonStore = 1000.0 	
    terrestrialBiosphereCarbonStore = 3000.0   
    atmosCarbonStore = [0 for i in range(years+1)]
    atmosCarbonStore[0]=atmosCarbonStoreStart

    # Flux coefficients
    atmosCarbonCoeff=0.15
    oceanDegassingCoeff = 1e-25
    oceanDegassingCoeffPow= 9.0    
    photosynthesisCoeff = 16.4 
    photosynthesisCoeffPow =0.2   
    respirationCoeff = 0.019                  
         
    # Initial values of fluxes (Gt/yr)
    oceanCarbonUptakeFlux = atmosCarbonStore[0] *atmosCarbonCoeff                 
    oceanDegassingFlux =  oceanDegassingCoeff * surfaceOceanCarbonStore**oceanDegassingCoeffPow
    photosynthesisFlux = photosynthesisCoeff  * atmosCarbonStore[0]**photosynthesisCoeffPow 
    respirationFlux = respirationCoeff * terrestrialBiosphereCarbonStore                   

    for i in range(years):

        atmosCarbonStore[i+1] = atmosCarbonStore[i] + (respirationFlux + oceanDegassingFlux + total_anthropogenic_emissions[i] -  photosynthesisFlux - oceanCarbonUptakeFlux) * dt
        surfaceOceanCarbonStore = surfaceOceanCarbonStore + (oceanCarbonUptakeFlux - oceanDegassingFlux) * dt
        terrestrialBiosphereCarbonStore = terrestrialBiosphereCarbonStore + (photosynthesisFlux - respirationFlux) * dt
        oceanCarbonUptakeFlux = atmosCarbonStore[i] * atmosCarbonCoeff             
        oceanDegassingFlux = oceanDegassingCoeff * surfaceOceanCarbonStore**oceanDegassingCoeffPow          
        photosynthesisFlux = photosynthesisCoeff * atmosCarbonStore[i]**photosynthesisCoeffPow           
        respirationFlux = respirationCoeff *terrestrialBiosphereCarbonStore                
 
     #Convert back to PPM.
    for i in range(years+1):
        atmosCarbonStore[i]=atmosCarbonStore[i]/PPM2GT

    return atmosCarbonStore

def plot_carbon(co2):

# CO2 concentration in ppm between 1959 and 2014
# from NOAA Mauna Loa observatory record from:
# http://co2now.org/Current-CO2/CO2-Now/noaa-mauna-loa-co2-data.html */
    co2AtmosConcentrationInPPM =[
    315.97,316.91,317.64,318.45,318.99,319.62,320.04,321.38,322.16,323.04,
    324.62,325.68,326.32,327.45,329.68,330.18,331.08,332.05,333.78,335.41,
    336.78,338.68,340.10,341.44,343.03,344.58,346.04,347.39,349.16,351.56,
    353.07,354.35,355.57,356.38,357.07,358.82,360.80,362.59,363.71,366.65,
    368.33,369.52,371.13,373.22,375.77,377.49,379.80,381.90,383.76,385.59,
    387.37,389.85,391.63,393.82,396.48,398.55]

    y=[]
    for i in range(1959,2015):
        y.append(i)

    xmarks=[i for i in range(1959,2015+1,5)]

    plt.xlabel("Year") 
    plt.ylabel("Atmospheric carbon dioxide concentrations (PPM)") 
    plt.xticks(xmarks,rotation=90) 
    plt.plot(y,co2, label="Model")
    plt.plot(y,co2AtmosConcentrationInPPM,label="Observations")
    
  
    plt.legend()

    fig = plt.gcf() 
    fig.subplots_adjust(bottom=0.2)
    fig.savefig("carbon.pdf", format='pdf',dpi=720)
    plt.close()

def plot_future():
    
    initialAtmosCarbon=330.35
    years=110
    scenarios=["A1", "A1T", "A1G", "A2", "B1", "B2"]

    y=[]
    for i in range(1990,2101):
        y.append(i)

    xmarks=[i for i in range(1990,2100+1,10)]

    plt.xlabel("Year") 
    plt.ylabel("Atmospheric carbon dioxide concentrations (PPM)") 
    plt.xticks(xmarks,rotation=90) 
    
    for i in scenarios:
        emissions=emissionsScenario(i)
        c=box_model(initialAtmosCarbon,years,emissions)
        plt.plot(y,c, label=i)

    plt.legend()
    fig = plt.gcf() 
    fig.subplots_adjust(bottom=0.2)
    #fig.set_size_inches(10, 8)
    fig.savefig("carbon_predict.pdf", format='pdf',dpi=720)
    plt.close()

def plot_emission_scenarios():
  
    years=110
    scenarios=["A1", "A1T", "A1G", "A2", "B1", "B2"]

    y=[]
    for i in range(1990,2100):
        y.append(i)

    xmarks=[i for i in range(1990,2100+1,10)]

    plt.xlabel("Year") 
    plt.ylabel("Future human emission scenarios of carbon (GT)") 
    plt.xticks(xmarks,rotation=90) 
    
    for i in scenarios:
        emissions=emissionsScenario(i)
        plt.plot(y,emissions, label=i)

    plt.legend()
    fig = plt.gcf() 
    fig.subplots_adjust(bottom=0.2)
    fig.savefig("carbon_emission_scenarios.pdf", format='pdf',dpi=720)
    plt.close()

def plot_emission_actual():

    emissions=calc_total_anthropogenic_emissions()

    y=[]
    for i in range(1959,2018):
        y.append(i)

    xmarks=[i for i in range(1959,2015+1,5)]

    plt.xlabel("Year") 
    plt.ylabel("Actual human emissions of carbon (GT)") 
    plt.xticks(xmarks,rotation=90)
    emissions=calc_total_anthropogenic_emissions()
    plt.plot(y,emissions, label="Total Emissions (GT)")
    emissions=calc_total_anthropogenic_emissions("FF")
    plt.plot(y,emissions, label="Fossil Fuel Emissions (GT)")
    emissions=calc_total_anthropogenic_emissions("LUC")
    plt.plot(y,emissions, label="Land Use Change Emissions (GT)")
    plt.legend()
    fig = plt.gcf() 
    fig.subplots_adjust(bottom=0.2)
    fig.savefig("carbon_emission_actual.pdf", format='pdf',dpi=720)
    plt.close()


initialAtmosCarbon=315.97 # in ppm
years=55
total_anthropogenic_emissions=calc_total_anthropogenic_emissions()
co2=box_model(initialAtmosCarbon,years,total_anthropogenic_emissions)

plot_carbon(co2)
plot_future()
plot_emission_scenarios()
plot_emission_actual()
