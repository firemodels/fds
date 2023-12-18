
# This script reads the experimental heat release rate and ventilation rates and appends RAMPs to the FDS input files

HRR_0 = 10000.0
T = [0.0] * 1000
V_DOT_NORTH = [0.0] * 1000
V_DOT_SOUTH = [0.0] * 1000

with open('paramfile.csv', 'r') as file:
    next(file)
    for line in file:
        INPUT_FILE, TEST_ID, T_END, TMPA, NOM_HRR, SEQ, N_HEIGHT, S_HEIGHT = line.strip().split(',')

        HRR_FILE = '../../../../exp/Memorial_Tunnel/HR' + TEST_ID.strip() + '.csv'
        with open(HRR_FILE, 'r') as hrr_file:
            next(hrr_file)
            next(hrr_file)
            with open(INPUT_FILE, 'a') as output_file:
                while True:
                    line = hrr_file.readline()
                    if not line:
                        output_file.write("\n")
                        break
                    TT, HRR, HRR2 = line.strip().split(',')
                    output_file.write("&RAMP ID='RAMP_FIRE', T={0:.1f}, F={1:.3f} /\n".format(float(TT), max(0.0, float(HRR) / HRR_0)))

        if int(SEQ) == 1:
            N_SOUTH_VENTS = 100
            N_NORTH_VENTS = 156
        elif int(SEQ) == 3:
            N_SOUTH_VENTS = 78
            N_NORTH_VENTS = 78
        elif int(SEQ) == 4:
            N_SOUTH_VENTS = 100
            N_NORTH_VENTS = 100
        elif int(SEQ) == 5:
            N_SOUTH_VENTS = 78
            N_NORTH_VENTS = 78
        elif int(SEQ) == 6:
            N_SOUTH_VENTS = 78
            N_NORTH_VENTS = 78
        elif int(SEQ) == 8:
            N_SOUTH_VENTS = 78
            N_NORTH_VENTS = 78
        elif int(SEQ) == 9:
            N_SOUTH_VENTS = 1
            N_NORTH_VENTS = 2
        elif int(SEQ) == 101:
            N_SOUTH_VENTS = 1
            N_NORTH_VENTS = 1
        elif int(SEQ) == 13:
            N_SOUTH_VENTS = 78
            N_NORTH_VENTS = 78
        elif int(SEQ) == 14:
            N_SOUTH_VENTS = 78
            N_NORTH_VENTS = 78

        VENT_FILE = '../../../../exp/Memorial_Tunnel/AP' + TEST_ID.strip() + '.csv'
        i=0
        with open(VENT_FILE, 'r') as vent_file:
            next(vent_file)
            next(vent_file)
            with open(INPUT_FILE, 'a') as output_file:
                while True:
                    line = vent_file.readline()
                    if not line:
                        output_file.write("\n")
                        break
                    i = i+1
                    T[i], V_DOT_SOUTH[i],V_DOT_NORTH[i] = line.strip().split(',')
                    output_file.write("&RAMP ID='RAMP_SOUTH', T={0:.1f}, F={1:.3f} /\n".format(float(T[i]), abs(float(V_DOT_SOUTH[i]) / float(N_SOUTH_VENTS))))
                for ii in range(1,i+1):
                    output_file.write("&RAMP ID='RAMP_NORTH', T={0:.1f}, F={1:.3f} /\n".format(float(T[ii]), abs(float(V_DOT_NORTH[ii]) / float(N_NORTH_VENTS))))

                output_file.write("\n")
                output_file.write("&TAIL /\n")
