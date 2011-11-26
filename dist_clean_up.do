*Clean up distance data
clear

cd "H://My Documents/Lucas-Alvarez-Replication"
use dist_cepii.dta

keep if (iso_o == "USA" | iso_o == "USA" | iso_o == "JPN" | iso_o == "JPN" | iso_o == "DEU" | iso_o == "DEU" | iso_o == "FRA" | iso_o == "FRA" | iso_o == "GBR" | iso_o == "GBR" | iso_o == "ITA" | iso_o == "ITA" | iso_o == "CHN" | iso_o == "CHN" | iso_o == "BRA" | iso_o == "BRA" | iso_o == "CAN" | iso_o == "CAN" | iso_o == "ESP" | iso_o == "ESP" | iso_o == "MEX" | iso_o == "MEX" | iso_o == "IND" | iso_o == "IND" | iso_o == "AUS" | iso_o == "AUS" | iso_o == "NLD" | iso_o == "NLD" | iso_o == "RUS" | iso_o == "RUS" | iso_o == "ARG" | iso_o == "ARG" | iso_o == "CHE" | iso_o == "CHE" | iso_o == "BEL" | iso_o == "BEL" | iso_o == "SWE" | iso_o == "SWE" | iso_o == "AUT" | iso_o == "AUT" | iso_o == "TUR" | iso_o == "TUR" | iso_o == "IDN" | iso_o == "IDN" | iso_o == "DNK" | iso_o == "DNK" | iso_o == "HKG" | iso_o == "HKG" | iso_o == "NOR" | iso_o == "NOR" | iso_o == "THA" | iso_o == "THA" | iso_o == "POL" | iso_o == "POL" | iso_o == "SAU" | iso_o == "SAU" | iso_o == "ZAF" | iso_o == "ZAF" | iso_o == "FIN" | iso_o == "FIN" | iso_o == "GRC" | iso_o == "GRC" | iso_o == "PRT" | iso_o == "PRT" | iso_o == "ISR" | iso_o == "ISR" | iso_o == "IRN" | iso_o == "IRN" | iso_o == "COL" | iso_o == "COL" | iso_o == "VEN" | iso_o == "VEN" | iso_o == "MYS" | iso_o == "MYS" | iso_o == "SGP" | iso_o == "SGP" | iso_o == "IRL" | iso_o == "IRL" | iso_o == "EGY" | iso_o == "EGY" | iso_o == "PHL" | iso_o == "PHL" | iso_o == "CHL" | iso_o == "CHL" | iso_o == "PAK" | iso_o == "PAK" | iso_o == "NZL" | iso_o == "NZL" | iso_o == "PER" | iso_o == "PER" | iso_o == "CZE" | iso_o == "CZE" | iso_o == "DZA" | iso_o == "DZA" | iso_o == "HUN" | iso_o == "HUN" | iso_o == "UKR" | iso_o == "UKR" | iso_o == "BGD" | iso_o == "BGD" | iso_o == "ROM" | iso_o == "ROM" | iso_o == "MAR" | iso_o == "MAR" | iso_o == "NGA" | iso_o == "NGA" | iso_o == "VNM" | iso_o == "VNM" | iso_o == "BLR" | iso_o == "BLR" | iso_o == "KAZ" | iso_o == "KAZ" | iso_o == "SVK" | iso_o == "SVK" | iso_o == "TUN" | iso_o == "TUN" | iso_o == "LKA" | iso_o == "LKA") & (iso_d == "USA" | iso_d == "USA" | iso_d == "JPN" | iso_d == "JPN" | iso_d == "DEU" | iso_d == "DEU" | iso_d == "FRA" | iso_d == "FRA" | iso_d == "GBR" | iso_d == "GBR" | iso_d == "ITA" | iso_d == "ITA" | iso_d == "CHN" | iso_d == "CHN" | iso_d == "BRA" | iso_d == "BRA" | iso_d == "CAN" | iso_d == "CAN" | iso_d == "ESP" | iso_d == "ESP" | iso_d == "MEX" | iso_d == "MEX" | iso_d == "IND" | iso_d == "IND" | iso_d == "AUS" | iso_d == "AUS" | iso_d == "NLD" | iso_d == "NLD" | iso_d == "RUS" | iso_d == "RUS" | iso_d == "ARG" | iso_d == "ARG" | iso_d == "CHE" | iso_d == "CHE" | iso_d == "BEL" | iso_d == "BEL" | iso_d == "SWE" | iso_d == "SWE" | iso_d == "AUT" | iso_d == "AUT" | iso_d == "TUR" | iso_d == "TUR" | iso_d == "IDN" | iso_d == "IDN" | iso_d == "DNK" | iso_d == "DNK" | iso_d == "HKG" | iso_d == "HKG" | iso_d == "NOR" | iso_d == "NOR" | iso_d == "THA" | iso_d == "THA" | iso_d == "POL" | iso_d == "POL" | iso_d == "SAU" | iso_d == "SAU" | iso_d == "ZAF" | iso_d == "ZAF" | iso_d == "FIN" | iso_d == "FIN" | iso_d == "GRC" | iso_d == "GRC" | iso_d == "PRT" | iso_d == "PRT" | iso_d == "ISR" | iso_d == "ISR" | iso_d == "IRN" | iso_d == "IRN" | iso_d == "COL" | iso_d == "COL" | iso_d == "VEN" | iso_d == "VEN" | iso_d == "MYS" | iso_d == "MYS" | iso_d == "SGP" | iso_d == "SGP" | iso_d == "IRL" | iso_d == "IRL" | iso_d == "EGY" | iso_d == "EGY" | iso_d == "PHL" | iso_d == "PHL" | iso_d == "CHL" | iso_d == "CHL" | iso_d == "PAK" | iso_d == "PAK" | iso_d == "NZL" | iso_d == "NZL" | iso_d == "PER" | iso_d == "PER" | iso_d == "CZE" | iso_d == "CZE" | iso_d == "DZA" | iso_d == "DZA" | iso_d == "HUN" | iso_d == "HUN" | iso_d == "UKR" | iso_d == "UKR" | iso_d == "BGD" | iso_d == "BGD" | iso_d == "ROM" | iso_d == "ROM" | iso_d == "MAR" | iso_d == "MAR" | iso_d == "NGA" | iso_d == "NGA" | iso_d == "VNM" | iso_d == "VNM" | iso_d == "BLR" | iso_d == "BLR" | iso_d == "KAZ" | iso_d == "KAZ" | iso_d == "SVK" | iso_d == "SVK" | iso_d == "TUN" | iso_d == "TUN" | iso_d == "LKA" | iso_d == "LKA")

gen con_num_o = 0
gen con_num_d = 0

replace con_num_o = 1 if iso_o == "USA"
replace con_num_o = 2 if iso_o == "JPN"
replace con_num_o = 3 if iso_o == "DEU"
replace con_num_o = 4 if iso_o == "FRA"
replace con_num_o = 5 if iso_o == "GBR"
replace con_num_o = 6 if iso_o == "ITA"
replace con_num_o = 7 if iso_o == "CHN"
replace con_num_o = 8 if iso_o == "BRA"
replace con_num_o = 9 if iso_o == "CAN"
replace con_num_o = 10 if iso_o == "ESP"
replace con_num_o = 11 if iso_o == "MEX"
replace con_num_o = 12 if iso_o == "IND"
replace con_num_o = 13 if iso_o == "AUS"
replace con_num_o = 14 if iso_o == "NLD"
replace con_num_o = 15 if iso_o == "RUS"
replace con_num_o = 16 if iso_o == "ARG"
replace con_num_o = 17 if iso_o == "CHE"
replace con_num_o = 18 if iso_o == "BEL"
replace con_num_o = 19 if iso_o == "SWE"
replace con_num_o = 20 if iso_o == "AUT"
replace con_num_o = 21 if iso_o == "TUR"
replace con_num_o = 22 if iso_o == "IDN"
replace con_num_o = 23 if iso_o == "DNK"
replace con_num_o = 24 if iso_o == "HKG"
replace con_num_o = 25 if iso_o == "NOR"
replace con_num_o = 26 if iso_o == "THA"
replace con_num_o = 27 if iso_o == "POL"
replace con_num_o = 28 if iso_o == "SAU"
replace con_num_o = 29 if iso_o == "ZAF"
replace con_num_o = 30 if iso_o == "FIN"
replace con_num_o = 31 if iso_o == "GRC"
replace con_num_o = 32 if iso_o == "PRT"
replace con_num_o = 33 if iso_o == "ISR"
replace con_num_o = 34 if iso_o == "IRN"
replace con_num_o = 35 if iso_o == "COL"
replace con_num_o = 36 if iso_o == "VEN"
replace con_num_o = 37 if iso_o == "MYS"
replace con_num_o = 38 if iso_o == "SGP"
replace con_num_o = 39 if iso_o == "IRL"
replace con_num_o = 40 if iso_o == "EGY"
replace con_num_o = 41 if iso_o == "PHL"
replace con_num_o = 42 if iso_o == "CHL"
replace con_num_o = 43 if iso_o == "PAK"
replace con_num_o = 44 if iso_o == "NZL"
replace con_num_o = 45 if iso_o == "PER"
replace con_num_o = 46 if iso_o == "CZE"
replace con_num_o = 47 if iso_o == "DZA"
replace con_num_o = 48 if iso_o == "HUN"
replace con_num_o = 49 if iso_o == "UKR"
replace con_num_o = 50 if iso_o == "BGD"
replace con_num_o = 51 if iso_o == "ROM"
replace con_num_o = 52 if iso_o == "MAR"
replace con_num_o = 53 if iso_o == "NGA"
replace con_num_o = 54 if iso_o == "VNM"
replace con_num_o = 55 if iso_o == "BLR"
replace con_num_o = 56 if iso_o == "KAZ"
replace con_num_o = 57 if iso_o == "SVK"
replace con_num_o = 58 if iso_o == "TUN"
replace con_num_o = 59 if iso_o == "LKA"

replace con_num_d = 1 if iso_d == "USA"
replace con_num_d = 2 if iso_d == "JPN"
replace con_num_d = 3 if iso_d == "DEU"
replace con_num_d = 4 if iso_d == "FRA"
replace con_num_d = 5 if iso_d == "GBR"
replace con_num_d = 6 if iso_d == "ITA"
replace con_num_d = 7 if iso_d == "CHN"
replace con_num_d = 8 if iso_d == "BRA"
replace con_num_d = 9 if iso_d == "CAN"
replace con_num_d = 10 if iso_d == "ESP"
replace con_num_d = 11 if iso_d == "MEX"
replace con_num_d = 12 if iso_d == "IND"
replace con_num_d = 13 if iso_d == "AUS"
replace con_num_d = 14 if iso_d == "NLD"
replace con_num_d = 15 if iso_d == "RUS"
replace con_num_d = 16 if iso_d == "ARG"
replace con_num_d = 17 if iso_d == "CHE"
replace con_num_d = 18 if iso_d == "BEL"
replace con_num_d = 19 if iso_d == "SWE"
replace con_num_d = 20 if iso_d == "AUT"
replace con_num_d = 21 if iso_d == "TUR"
replace con_num_d = 22 if iso_d == "IDN"
replace con_num_d = 23 if iso_d == "DNK"
replace con_num_d = 24 if iso_d == "HKG"
replace con_num_d = 25 if iso_d == "NOR"
replace con_num_d = 26 if iso_d == "THA"
replace con_num_d = 27 if iso_d == "POL"
replace con_num_d = 28 if iso_d == "SAU"
replace con_num_d = 29 if iso_d == "ZAF"
replace con_num_d = 30 if iso_d == "FIN"
replace con_num_d = 31 if iso_d == "GRC"
replace con_num_d = 32 if iso_d == "PRT"
replace con_num_d = 33 if iso_d == "ISR"
replace con_num_d = 34 if iso_d == "IRN"
replace con_num_d = 35 if iso_d == "COL"
replace con_num_d = 36 if iso_d == "VEN"
replace con_num_d = 37 if iso_d == "MYS"
replace con_num_d = 38 if iso_d == "SGP"
replace con_num_d = 39 if iso_d == "IRL"
replace con_num_d = 40 if iso_d == "EGY"
replace con_num_d = 41 if iso_d == "PHL"
replace con_num_d = 42 if iso_d == "CHL"
replace con_num_d = 43 if iso_d == "PAK"
replace con_num_d = 44 if iso_d == "NZL"
replace con_num_d = 45 if iso_d == "PER"
replace con_num_d = 46 if iso_d == "CZE"
replace con_num_d = 47 if iso_d == "DZA"
replace con_num_d = 48 if iso_d == "HUN"
replace con_num_d = 49 if iso_d == "UKR"
replace con_num_d = 50 if iso_d == "BGD"
replace con_num_d = 51 if iso_d == "ROM"
replace con_num_d = 52 if iso_d == "MAR"
replace con_num_d = 53 if iso_d == "NGA"
replace con_num_d = 54 if iso_d == "VNM"
replace con_num_d = 55 if iso_d == "BLR"
replace con_num_d = 56 if iso_d == "KAZ"
replace con_num_d = 57 if iso_d == "SVK"
replace con_num_d = 58 if iso_d == "TUN"
replace con_num_d = 59 if iso_d == "LKA"

keep con_num_o con_num_d distcap
outfile using dist_dat.txt, replace c

