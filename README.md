# radiosondes
tool for calculating pressure from temperature, relative humidity and GPS data

### Usage
you can install this package in a new environment or use the parent-folder (the one containing setup.py) as working directory.
I would rather use the working directory approach so far.

`from datetime import datetime`\
`from radiosondes.radiosondes.radiosondes import RS_Trainou2019`\
`path = r"<path-to-your-radiosonde-input-file"`\
`P0 = <launch pressure>`\
`z0 = <launch altitude>`\
`mirror = False`\
`P_Calc = True`\
`ipol = True`\
`geometric = False` #don't use the original altitude, but the calculated geopotential altitude\
`dirksenTv = True` #this calculation is faster than if set to `True`\
`start = datetime(<year>,<month>,<day>,<hour>,<min>,<sec>)`\

`RS = RS_Trainou2019(path=path, delim=' ', P0=P0, z0=z0, mirror=mirror, P_Calc=P_Calc, ipol=ipol, geometric=False, dirksenTv=True, start=start)`

Now you created a RS_Trainou2019 object with an updated dataframe
