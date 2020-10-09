EESchema Schematic File Version 4
LIBS:tri-slot-feedlines-cache
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
$Comp
L Device:C C1
U 1 1 5DF34BA3
P 5500 3050
F 0 "C1" V 5248 3050 50  0000 C CNN
F 1 "C" V 5339 3050 50  0000 C CNN
F 2 "Capacitor_SMD:C_0805_2012Metric_Pad1.15x1.40mm_HandSolder" H 5538 2900 50  0001 C CNN
F 3 "~" H 5500 3050 50  0001 C CNN
	1    5500 3050
	0    1    1    0   
$EndComp
$Comp
L Device:C C2
U 1 1 5DF34C46
P 5500 3600
F 0 "C2" V 5248 3600 50  0000 C CNN
F 1 "C" V 5339 3600 50  0000 C CNN
F 2 "Capacitor_SMD:C_0805_2012Metric_Pad1.15x1.40mm_HandSolder" H 5538 3450 50  0001 C CNN
F 3 "~" H 5500 3600 50  0001 C CNN
	1    5500 3600
	0    1    1    0   
$EndComp
$Comp
L Connector:TestPoint slot_l1
U 1 1 5DF34D0F
P 3950 3050
F 0 "slot_l1" H 4008 3170 50  0000 L CNN
F 1 "TestPoint" H 3500 3150 50  0001 L CNN
F 2 "TestPoint:TestPoint_Loop_D3.80mm_Drill2.0mm" H 4150 3050 50  0001 C CNN
F 3 "~" H 4150 3050 50  0001 C CNN
	1    3950 3050
	1    0    0    -1  
$EndComp
Wire Wire Line
	3950 3050 4800 3050
Wire Wire Line
	3950 3600 4800 3600
Wire Wire Line
	5650 3050 7050 3050
Wire Wire Line
	5650 3600 7050 3600
Text Label 7050 3050 0    50   ~ 0
feed_l
Text Label 7050 3600 0    50   ~ 0
feed_r
$Comp
L Connector:TestPoint slot_r1
U 1 1 5DF351B5
P 3950 3600
F 0 "slot_r1" H 4008 3720 50  0000 L CNN
F 1 "TestPoint" H 3500 3700 50  0001 L CNN
F 2 "TestPoint:TestPoint_Loop_D3.80mm_Drill2.0mm" H 4150 3600 50  0001 C CNN
F 3 "~" H 4150 3600 50  0001 C CNN
	1    3950 3600
	1    0    0    -1  
$EndComp
$Comp
L Device:C C3
U 1 1 5DF4506D
P 4800 3300
F 0 "C3" V 4548 3300 50  0000 C CNN
F 1 "C" V 4639 3300 50  0000 C CNN
F 2 "Capacitor_SMD:C_0805_2012Metric_Pad1.15x1.40mm_HandSolder" H 4838 3150 50  0001 C CNN
F 3 "~" H 4800 3300 50  0001 C CNN
	1    4800 3300
	-1   0    0    1   
$EndComp
Wire Wire Line
	4800 3150 4800 3050
Connection ~ 4800 3050
Wire Wire Line
	4800 3050 5350 3050
Wire Wire Line
	4800 3450 4800 3600
Connection ~ 4800 3600
Wire Wire Line
	4800 3600 5350 3600
$EndSCHEMATC
