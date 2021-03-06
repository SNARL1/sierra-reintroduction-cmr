# Data from: Disease and climate effects on individuals drive post-reintroduction population dynamics of an endangered amphibian

*Maxwell B. Joseph and Roland A. Knapp*

This repository contains data to understand how the amphibian 
chytrid fungus *Batrachochytrium dendrobatidis* and climatic conditions drive 
population dynamics of Sierra Nevada yellow-legged frogs.
It combines a capture mark-recapture model with a joint model for Bd/host dynamics that accomodates known introductions of hosts.
This is associated with the paper **Disease and climate effects on individuals jointly drive post-reintroduction population dynamics of an endangered amphibian**.
This repo is maintained by Max Joseph (maxwell.b.joseph (at) colorado.edu), 
but questions about the data should be directed to both Max Joseph and 
Roland Knapp (roland.knapp (at) ucsb.edu).

The code associated with this project are available on GitHub: 
https://www.github.com/snarl1/sierra-mr

All of the raw data for this project are tabular, in comma separated value 
(CSV) format.
Missing values for all files are left blank (mark-recapture data) or coded as 
NA (snow data). 
There are seven CSV files:

```
.
+-- dana_meadows_snow_depth_2006-2017.csv
+-- dana_meadows_swe_2006-2017.csv
+-- ellery_lake_temperature_2006-2017.csv
+-- alpine_captures_2006-2017.csv
+-- subalpine_captures_2006-2017.csv
+-- alpine_introductions_2006-2017.csv
+-- alpine_introductions_2006-2017.csv
```

### Weather data

#### dana_meadows_snow_depth_2006-2017.csv

This is a csv file with daily snow depth as measured at Dana Meadows, which is 
used to model survey occurrence.
It was generated using data acquired from the California Data Exchange: http://cdec.water.ca.gov/cgi-progs/getDailyCSV?station_id=DAN&dur_code=D&sensor_num=18&start_date=2006/01/01&end_date=2017/12/31

Fields: 

- `date`: date in YYYY-MM-DD format
- `snow_depth_inches`: daily snow depth in inches

#### dana_meadows_swe_2006-2017.csv

This is a csv file with yearly snow water equivalent (SWE) as measured at Dana 
Meadows at the end of March or beginning of April (depending on year).
It was generated using data acquired from the California Data Exchange: http://cdec.water.ca.gov/cgi-progs/snowQuery?course_num=DAN&month=April&start_date=2006&end_date=2017&csv_mode=Y&data_wish=Retrieve+Data

Fields: 

- `measure_date`: date that SWE was measured in YYYY-MM-DD format
- `year`: the year 
- `snow_water_equivalent`: snow water equivalent for the year

#### ellery_lake_temperature_2006-2017.csv

This csv file contains hourly air temperature data measured at Ellery Lake, 
measured between 0900 and 1600 (not inclusive). 
This was generated by the California Data Exchange: 
http://cdec.water.ca.gov/cgi-progs/queryCSV?station_id=ery&sensor_num=4&dur_code=H&start_date=2006-06-01&end_date=2017-09-30&data_wish=View+CSV+Data

Fields: 

- `time`: time of day (military time format)
- `date`: measurement date (YYYY-MM-DD)
- `temp_c`: air temperature in degrees C

### Mark-recapture data

#### alpine_captures_2006-2017.csv and subalpine_captures_2008-2017.csv

These files contain the capture/recapture data for the Alpine and Subalpine 
study sites. 

Fields: 

- `site_id`: Yosemite National Park site identification code
- `capture_date`: date of capture YYYY-MM-DD
- `survey_year`: year of survey
- `period_id`: survey period id, in format `a.b.c`, where `a` is a numeric index for study year, starting at 1, `b` is a primary period index within years, and `c` is a secondary period or survey index within primary periods
- `pit_tag_id`: unique identifier associated with an individual's PIT tag 
- `swab_id`: unique identification code associated with each skin swab sample
- `bd_load`: number of *Batrachochytrium dendrobatidis* (Bd) ITS copies on the swab

#### alpine_introductions_2006-2017.csv and subalpine_introductions_2008-2017.csv

These files contain information on which individuals were translocated 
(introduced) into the Alpine and Subalpine populations, and when. 

Fields: 

- `site_id`: Yosemite National Park site identification code
- `translocate_date`: the date on which an individual was introduced YYYY-MM-DD
- `pit_tag_id`: unique identifier associated with an individual's PIT tag 

\* A note on study site pseudonyms: We have anonymized the study sites because 
they harbor endangered species. 
But, we have retained the numeric site codes in the data. 
These can be back-referenced to actual locations by Yosemite National Park. 
