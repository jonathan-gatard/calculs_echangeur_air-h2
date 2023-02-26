Project ISAE SUPAERO/ENSMA 2018

This code performs calculations for a heat exchanger that consists of tubes with a certain diameter and length, between dihydrogen and air.
First I initialize various parameters such as the dimensions of the heat exchanger, temperatures, flow rates, and material properties.
Then I proceed to calculate the number of tubes that can fit on each level of the heat exchanger and the number of levels that can fit in the heat exchanger based on the given parameters.

I use several functions to calculate various properties of air and dihydrogen such as viscosity, thermal conductivity, and specific heat.
These properties are used to calculate the Reynolds number and Nusselt number for each fluid in each tube.
The heat transfer coefficient is then calculated based on the Nusselt number and used to calculate the heat transfer rate for each tube.

The code also calculates the average temperatures of the air and dihydrogen at some point along the length of each tube.
The calculations are performed using iterative methods until a certain level of precision is achieved.

Finally, the code outputs the heat transfer rate for each tube.
