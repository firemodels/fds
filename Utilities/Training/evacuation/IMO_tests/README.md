# FDS+Evac | IMO Test Cases

This folder contains FDS+Evac input files of suggested tests that should be included in the verification process. The test cases defined by IMO MSC/Circ.1238 (30 October 2007) 'Guidelines For Evacuation Analyses for New and Existing Passenger Ships', Annex 3.

## Component testing

### Test 1: Maintaining set walking speed in bended_corridor

- IMO_Test_01.fds

### Test 2: Maintaining set walking speed up staircase

- IMO_Test_02_CORR.fds
- IMO_Test_02_EVSS.fds
- IMO_Test_02_EVSS2.fds

### Test 3: Maintaining set walking speed down staircase

### Test 4: Exit flow rate

Two different human parameter sets, `Adult 1` and `Adult 2`, are used and the resulting human flows are compared. Adult 1 is the standard 'Adult' defined in FDS+Evac, Adult 2 has `L_NON_SP=0.3`.

- IMO_Test_04a.fds
- IMO_Test_04b.fds

### Test 5: Response time

Response time, uniform 10-100 s.

- IMO_Test_05.fds

### Test 6: Rounding corners

- IMO_Test_06.fds

### Test 7: Assignment of population demographics parameters

The gender and age group dependent walking speeds given in the IMO. The demographics parameters of Annex 2 Table 3.4 are also tested.

- IMO_Test_07.fds
- IMO_Test_07_A2T3-4.fds

## Functional verification

## Qualitative verification

### Test 8: Counterflow - two rooms connected via a corridor

- IMO_Test_08a.fds
- IMO_Test_08b.fds
- IMO_Test_08c.fds
- IMO_Test_08d.fds

### Test 9: Exit flow: crowd dissipation from a large public room

- IMO_Test_09a.fds
- IMO_Test_09b.fds
- IMO_Test_09c.fds
- IMO_Test_09d.fds

### Test 10: Exit route allocation

- IMO_Test_10.fds

### Test 11: Staircase

- IMO_Test_11.fds

## Quantitative verification
