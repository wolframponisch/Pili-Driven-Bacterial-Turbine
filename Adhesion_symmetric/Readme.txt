%% To run the script, use below code:

# Set the parameters
TOTAL_TIME=10000
TIME_STEP=0.00001
CELL_NO=40
PILI_NO=14
RATE_T_ATT=2
TIME_T_DET=50
FORCE_T_DET=180
N_T_LENGTH=12
N_T_WIDTH=2
SEED=105

	g++ -O3 -o final.out main.cpp -std=c++0x -lm
	./final.out "$TOTAL_TIME" "$TIME_STEP" "$CELL_NO" "$PILI_NO" "$RATE_T_ATT" "$TIME_T_DET" "$FORCE_T_DET" "$N_T_LENGTH" "$N_T_WIDTH" "$SEED"

