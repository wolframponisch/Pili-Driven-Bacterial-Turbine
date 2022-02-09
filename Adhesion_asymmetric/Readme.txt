
# Set the parameters
        TOTAL_TIME=10000
        TIME_STEP=0.00001
        CELL_NO=40
        PILI_NO=15
        RATE_T_ATT=2        	# Attachment rate to turbine    
        RATE_T_ATT_SPECIAL=2	# Attachment rate to manipulated turbine region
        TIME_T_DET=2		# Detachment time from turbine
        TIME_T_DET_SPECIAL=50	# Detachment time from manipulated turbine region
        FORCE_T_DET=180
        N_T_LENGTH=12
        N_T_WIDTH=2
        SEED=281



g++ -O3 -o final.out main.cpp -std=c++0x -lm
./final.out "$TOTAL_TIME" "$TIME_STEP" "$CELL_NO" "$PILI_NO" "$RATE_T_ATT" "$RATE_T_ATT_SPECIAL" "$TIME_T_DET" "$TIME_T_DET_SPECIAL" "$FORCE_T_DET" "$N_T_LENGTH" "$N_T_WIDTH" "$SEED"
 

