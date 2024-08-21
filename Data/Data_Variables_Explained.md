# Data Variables Explained:

This is an accessory document to explain what each variable in the two datasets ("cumulative_temperatures_data" and "max_rates_data") means and/or how derived variables were calculated. For further information on how these variables were measures/calculated or their relevance to the analysis, please refer to the published paper (see full reference in this repository's *README* file).

## Cumulative Temperatures Dataset:

This dataset ("*cumulative_temperatures_data*") contains the cumulative temperature change (in ºC) of each body part (snout, eye, head, dorsum, leg, foot, tail) for each individual, in each trial of each treatment. Cumulative temperature change of each body part was calculated as: **temperature of body part *x* at time *y* - temperature of body part *x* at time=0**. These are then logged in the collumn headed by the name of the corresponding body part (column names: "Snout", "Eye", "Head", "Dorsum", "Leg", "Foot", "Tail"). 
The remaining variables on the dataset represent:

  - *"**treat**"*: The trial's treamtment as a 2 level categorical variable (Heliothermy vs Thigmothermy);

  - *"**heat**"*: The trial's direction of heat eaxchage (i.e. whether the animal was gaining or loosing heat) as a 2 level categorical variable (Heating vs Cooling);
    
  - *"**id**"*: The alphanumeric individual identification code given to each individual animal in the trial;

  - *"**time**"*: The trial's elapsed time, in seconds;

  - *"**t_water**"*: The temperature (in ºC) of the water in the experimental set-up, measured with a contact thermometer with k-type thermocouple probe;

  - *"**t_air**"*: The ambient (room) temperature (in ºC), measured with a temperature and humidity meter;

  - *"**posture**"*: The posture of the animal at the time of reading the remainder of the data (every 20 seconds), as a categorical variable with 3 mutually exclusive leves:
    - *0* -> the animal was laying completely flat against the surface of the arena
    - *0.5* -> half of the animal was flat against the arena and half was raised
    - *1* -> the animal's torsum and pelvis were fully raised (i.e. body raised up on all four legs) from the arena
    
    In all cases, the head and tail positions were disregarded.
    
  - *"**urine**"*: An account on whether the animal performed any excretions (water, urates and/or faeces) in the prior 20 seconds leading to that measurement;

  - *"**mass_initial**"*: The mass of the animal measured (with a precision balance) immediately **before** the start of that specific treatment (Heliothermy vs Thigmothermy) set of trials (i.e. heating trial and subsequent cooling trial);

  - *"**mass_initial**"*: The mass of the animal measured (with a precision balance) immediately **after** the end of that specific treatment (Heliothermy vs Thigmothermy) set of trials (i.e. heating trial and subsequent cooling trial);

  - *"**svl**"*: The Snout-to-Vent length of the animal measured with a digital calliper, before the start of the experiments;

  - *"**population**"*: The location/populaiton from which the animal was captured (and later returned to, after the end of the trials);

  - *"**feed**"*: The number of days since the animal was last fed (minimum was 1 day, since trials were always performed with animals on a post-absorvative state);

  - *"**jump**"*: Whether in the previous 20 second period, the animal made any attempt to escape the trial arena and/or the thermal camera's field-of-view. If yes then it was scored as *1*, if not the cell was left empty (*NA*);

  - *"**heat2**"*: A variant of the *"**heat**"* varible that simply assigns the correct order of heat exchange processes (animals first underwhent the heating trial and then, after that, the respective cooling trial) - i.e. 1=Heating, 2=Cooling. Simply derived to facilitate ordering the levels of the *"**heat**"* variable when plotting the results in *R*.     


## Maximum Rates Dataset:

This dataset ("*max_rates_data*") contains the cumulative temperature change (in ºC) of each body part (snout, eye, head, dorsum, leg, foot, tail) for each individual, in each trial of each treatment, **for a 1-minute period of between *time=60s* and *time=120s***. Cumulative temperature change of each body part was calculated as: **temperature of body part *x* at time *y* - temperature of body part *x* at time=60** and logged in the collumn headed *"**temp**"*. This period of time was deemed the most representative of the maximum rate of heat exchange due to the largest temperature difference between the animal and the environment. Time period 0-59s was disregarded since, being the start of the trial, animals often exhibited exploratory instead of thermoregulatory behaviours. 
Some variables (*"**treat**"*, *"**id**"* are the same as for the "*cumulative_temperature_data" dataset, while new ones represent:

  - *"**treat**"*: The
