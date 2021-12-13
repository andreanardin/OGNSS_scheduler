O-GNSS SCHEDULER
================

This file is part of OGNSS_scheduler.

OGNSS_scheduler is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OGNSS_scheduler is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with OGNSS_scheduler.  If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2020 Andrea Nardin
Navigation, Signal Analysis and Simulation (NavSAS) group,
Politecnico di Torino
andrea.nardin3@gmail.com


Description
-----------

This is a PDOP-driven Scheduler for Optical Inter-Satellite Links enabled GNSS.
A thorough description of the working principle can be found in:

A. Nardin, J. A. Fraire, F. Dovis, "Contact Plan Design for GNSS Constellations: A Case Study with Optical Inter-Satellite Links", IEEE Transactions on Aerospace and Electronic Systems

Please cite the relative paper if you use this software.

The scheduler is mainly composed by the following blocks:
1) Data retrieval from AGI STK access reports
2) Instantaneous PDOP computation for each potential schedulable link
3) Time slotting and PDOP averaging. Graph weights computation
4) Degree-constrained minimum spanning tree extraction
5) Completing available links assignation with respect to user selected criteria
6) Output data visualization and metrics evaluation


Usage
-----

Modify the settings as desired by acting on the file *initSettings.m*.
Launch the scheduler by running *OGNSS_scheduler.m*


Output Data and Metrics
-----------------------

### CSV output files

#### Data Delivery Latency
- *longestPaths.csv*
	contains the longest Meo-to-meo path for each time window
- *longestPathIds.csv*
	contains the indices of the MEO pair corresponding to the longst path, index ranges from 1 to Nmeo (MEO indices of the weight matrices are from Nleo+1 to Nleo+Nmeo instead)
- *averageDistances.csv*
	contains in each element (i,j) the average distance from MEO satellite i to any other MEO, at epoch j

#### Assignation Fairness
- *meoTotContactTime.csv*
	contains the total contact time windows assigned to each MEO satellite for the whole simulation. Each time window contact is accounted as 1, and lasted for tmin. No distinction is made between tuseful and tpoint, since the latter involves energy consumption as well.
- The total contact time is used to compute Jain Fairness Index included in table 1

#### Tables
- *table1.csv* contains:
	- **Disconnected ep.** the total number of epochs where the graph is not connected (should be zero)
	- **Avg idle antennas** The average number of idle antennas over the entire simulation (zero if secondary optimization is in place)
	- **No. of epochs with at least 1 B3MP** The number of epochs that have at least 1 LEO simultaneously connected to 3 different MEO orbital planes
	- **Max revisit time** The maximum revisit time experienced by a MEO during the simulation
	- **Ep. exceeding max node degree** The number of epochs where the maximum degree allowed for any node (i.e. the amount of onboard LCTs) is exceeded (should be zero)
	- **Longest distance** The longest distance (shortest path) ever experienced by any MEO pair over the simulation
	- **Longest distance ep.** The time window (epoch) when the longest distance (shortest path) is experienced.
	- **JFI** is the Jain Fairness Index, computed over th etotal contact time experienced by each MEO satellite.

- *table2.csv* contains:
	- **Percentiles** The selected PDOP percentiles that are evaluated
	- **PDOP** The PDOP values corresponing to each percentile



