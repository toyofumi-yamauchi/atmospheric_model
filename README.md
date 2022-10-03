# atmospheric_model
Import the data from NRLMSIS 2.0, calculate the atmospheric parameters

Author: Toyofumi Yamauchi
Last edit: October 03, 2022

If you want to make changes, please use a pull request.

How to use:
1.  Go to the webstite
    https://kauai.ccmc.gsfc.nasa.gov/instantrun/msis
2.  On the website, choose Input Parameters:
    Time Type = Universal
    Date + Time = 2000/01/01 01:30 (your choice)
    Coordinate Type = Geographic
    Latitude = 55 (your choice)
    Longitude = 45 (your choice)
    Height = 100 (your choice)
    Profile Type = Height [0 - 1000 km] (your choice)
    Start = 0 (your choice)
    Stop = 1000 (your choice)
    Step Size = 50 (your choice)
    F10.7 (daily) = -1
    F10.7 (3 month avg) = -1
    AP = -1
3.  On the website, choose Calculated MSIS Model Parameters
    Check everything (your choice)
    NRMSIS Version = NRLMSIS 2.0
    Output Option = Download Data
4.  Move the downloaded .txt file into this folder
5.  Open the python file
    - Mass_flow_rate.py to calculate the mass flow rate
    - Pressure.py to calculate the pressure
6.  Change the variable filename to the name you downloaded from the website
7.  Run the python file
    I use the Jupyter. 

