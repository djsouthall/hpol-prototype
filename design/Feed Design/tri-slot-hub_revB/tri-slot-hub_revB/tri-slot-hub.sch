EESchema Schematic File Version 4
LIBS:tri-slot-hub-cache
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
Text Label 7250 4300 0    50   ~ 0
slot_l
Text Label 6800 5550 0    50   ~ 0
slot_r
NoConn ~ 6450 6150
$Comp
L Connector:TestPoint TP1
U 1 1 5F2C7D1F
P 6250 4300
F 0 "TP1" H 6308 4420 50  0000 L CNN
F 1 "TestPoint" H 6308 4329 50  0000 L CNN
F 2 "TestPoint:TestPoint_Pad_1.0x1.0mm" H 6450 4300 50  0001 C CNN
F 3 "~" H 6450 4300 50  0001 C CNN
	1    6250 4300
	1    0    0    -1  
$EndComp
Wire Wire Line
	6250 4300 7250 4300
$Comp
L Connector:TestPoint TP2
U 1 1 5F2C7DAA
P 5700 4300
F 0 "TP2" H 5758 4420 50  0000 L CNN
F 1 "TestPoint" H 5758 4329 50  0000 L CNN
F 2 "TestPoint:TestPoint_Pad_1.0x1.0mm" H 5900 4300 50  0001 C CNN
F 3 "~" H 5900 4300 50  0001 C CNN
	1    5700 4300
	1    0    0    -1  
$EndComp
Wire Wire Line
	6250 4300 5700 4300
Connection ~ 6250 4300
$Comp
L Connector:TestPoint TP4
U 1 1 5F2C7E13
P 6250 4700
F 0 "TP4" H 6308 4820 50  0000 L CNN
F 1 "TestPoint" H 6308 4729 50  0000 L CNN
F 2 "TestPoint:TestPoint_Pad_1.0x1.0mm" H 6450 4700 50  0001 C CNN
F 3 "~" H 6450 4700 50  0001 C CNN
	1    6250 4700
	1    0    0    -1  
$EndComp
$Comp
L Connector:TestPoint TP3
U 1 1 5F2C7E1A
P 5700 4700
F 0 "TP3" H 5758 4820 50  0000 L CNN
F 1 "TestPoint" H 5758 4729 50  0000 L CNN
F 2 "TestPoint:TestPoint_Pad_1.0x1.0mm" H 5900 4700 50  0001 C CNN
F 3 "~" H 5900 4700 50  0001 C CNN
	1    5700 4700
	1    0    0    -1  
$EndComp
Wire Wire Line
	5700 4700 6250 4700
Connection ~ 6250 4700
Wire Wire Line
	6250 4700 6800 4700
$EndSCHEMATC
