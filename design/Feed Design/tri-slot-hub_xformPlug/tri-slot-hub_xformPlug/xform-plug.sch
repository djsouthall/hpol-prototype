EESchema Schematic File Version 4
LIBS:xform-plug-cache
EELAYER 26 0
EELAYER END
$Descr A4 11693 8268
encoding utf-8
Sheet 1 1
Title ""
Date ""
Rev ""
Comp ""
Comment1 ""
Comment2 ""
Comment3 ""
Comment4 ""
$EndDescr
Wire Wire Line
	6800 4700 6800 5550
Wire Wire Line
	6400 4300 7250 4300
Text Label 7250 4300 0    50   ~ 0
slot_l
Text Label 6800 5550 0    50   ~ 0
slot_r
$Comp
L Transformer:ADT4-6T TR1
U 1 1 5F2C57B8
P 6200 4500
F 0 "TR1" H 6200 4878 50  0000 C CNN
F 1 "ADT4-6T" H 6200 4787 50  0000 C CNN
F 2 "Package_SO:Mini-Circuits_CD637_H5.23mm" H 6200 4150 50  0001 C CNN
F 3 "https://ww2.minicircuits.com/pdfs/ADT4-6T+.pdf" H 6200 4500 50  0001 C CNN
	1    6200 4500
	1    0    0    -1  
$EndComp
Wire Wire Line
	6800 4700 6400 4700
$Comp
L eo_sch_symbols:GND #PWR0101
U 1 1 5F2C58BD
P 5300 4600
F 0 "#PWR0101" H 5300 4350 50  0001 C CNN
F 1 "GND" H 5305 4427 50  0000 C CNN
F 2 "" H 5300 4600 50  0001 C CNN
F 3 "" H 5300 4600 50  0001 C CNN
	1    5300 4600
	1    0    0    -1  
$EndComp
Wire Wire Line
	5300 4600 5300 4500
$Comp
L eo_Resistors:R0603 R2
U 1 1 5F2C594A
P 5900 5000
F 0 "R2" H 5968 5046 50  0000 L CNN
F 1 "0" H 5968 4955 50  0000 L CNN
F 2 "eo-footprints:0603_R" V 6250 5000 50  0001 C CNN
F 3 "~" V 6000 5000 50  0001 C CNN
F 4 "RES 1%  0603 (1608 Metric) " V 6400 5150 50  0001 C CNN "Description"
	1    5900 5000
	1    0    0    -1  
$EndComp
Wire Wire Line
	5900 4700 6000 4700
Wire Wire Line
	5900 4700 5900 4850
$Comp
L eo_sch_symbols:GND #PWR0102
U 1 1 5F2C5A17
P 5900 5250
F 0 "#PWR0102" H 5900 5000 50  0001 C CNN
F 1 "GND" H 5905 5077 50  0000 C CNN
F 2 "" H 5900 5250 50  0001 C CNN
F 3 "" H 5900 5250 50  0001 C CNN
	1    5900 5250
	1    0    0    -1  
$EndComp
Wire Wire Line
	5900 5250 5900 5150
$Comp
L eo_Resistors:R0603 R1
U 1 1 5F2C5A68
P 6800 4500
F 0 "R1" V 6900 4500 50  0000 C CNN
F 1 "0" V 6686 4500 50  0000 C CNN
F 2 "eo-footprints:0603_R" V 7150 4500 50  0001 C CNN
F 3 "~" V 6900 4500 50  0001 C CNN
F 4 "RES 1%  0603 (1608 Metric) " V 7300 4650 50  0001 C CNN "Description"
	1    6800 4500
	0    1    1    0   
$EndComp
Wire Wire Line
	6650 4500 6400 4500
$Comp
L eo_sch_symbols:GND #PWR0103
U 1 1 5F2C5B2A
P 7200 4500
F 0 "#PWR0103" H 7200 4250 50  0001 C CNN
F 1 "GND" V 7205 4372 50  0000 R CNN
F 2 "" H 7200 4500 50  0001 C CNN
F 3 "" H 7200 4500 50  0001 C CNN
	1    7200 4500
	0    -1   -1   0   
$EndComp
Wire Wire Line
	7200 4500 6950 4500
$Comp
L eo_connectors:CONSMA002-L J1
U 1 1 5F2C68F8
P 5300 4300
F 0 "J1" H 5376 4525 50  0000 C CNN
F 1 "CONSMA002-L" H 5376 4434 50  0000 C CNN
F 2 "eo-footprints:RF_SMA_RightAngle" H 5500 4500 60  0001 L CNN
F 3 "https://www.te.com/commerce/DocumentDelivery/DDEController?Action=srchrtrv&DocNm=1814400&DocType=Customer+Drawing&DocLang=English" H 5500 4600 60  0001 L CNN
F 4 "Connectors, Interconnects" H 5500 4900 60  0001 L CNN "Category"
F 5 "Coaxial Connectors (RF)" H 5500 5000 60  0001 L CNN "Family"
F 6 " CONN SMA RCPT STR 50OHM RA" H 5500 5300 60  0001 L CNN "Description"
F 7 "Active" H 5500 5500 60  0001 L CNN "Status"
F 8 "CONSMA003-L " H 5750 4850 50  0001 C CNN "Part Number"
F 9 "Linx Tech." H 5850 4950 50  0001 C CNN "Manufacturer"
	1    5300 4300
	1    0    0    -1  
$EndComp
Wire Wire Line
	5500 4300 6000 4300
$EndSCHEMATC
