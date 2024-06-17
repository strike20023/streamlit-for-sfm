#!/bin/bash

TIMESTEPS=1000000
EPISODES=10
RESET_NUM_TIMESTEPS=False
ENT_COEF=0.01
LEARNING_RATE=0.0003
N_STEPS=2048
BATCH_SIZE=64
N_EPOCHS=10
GAMMA=0.99
CLIP_RANGE=0.20
MAX_GRAD_NORM=0.5
VERBOSE=1
RANDOM_SEED=42
DEVICE="cuda"
PPO_POLICY="MlpPolicy"
REWARD_THRESHOLD=0.05
OUTPUT_FILE_PATH="Results/state_reward_6-5_5.31.txt"
RESULTS_PATH="Results/results_6-5_5.31.csv"
DATA_DIR="/mnt/ST8000/jialei/SFM_E/data/6-5"
VARIANCE=0.1
VARIANCE_MIN=0.3
VARIANCE_MAX=2.1
PERCENTAGE=0.05
PERCENTAGE_MIN=0.0
PERCENTAGE_MAX=1.0
FXY_ACTION=0.02
FXY_ACTION_MIN=0.0
FXY_ACTION_MAX=1.0
CXY_ACTION=0.05
CXY_ACTION_MIN=0.0
CXY_ACTION_MAX=1.0

python main.py \
  --TIMESTEPS "$TIMESTEPS" \
  --episodes "$EPISODES" \
  --reset_num_timesteps "$RESET_NUM_TIMESTEPS" \
  --ent_coef "$ENT_COEF" \
  --learning_rate "$LEARNING_RATE" \
  --n_steps "$N_STEPS" \
  --batch_size "$BATCH_SIZE" \
  --n_epochs "$N_EPOCHS" \
  --gamma "$GAMMA" \
  --clip_range "$CLIP_RANGE" \
  --max_grad_norm "$MAX_GRAD_NORM" \
  --verbose "$VERBOSE" \
  --random_seed "$RANDOM_SEED" \
  --device "$DEVICE" \
  --ppo_Policy "$PPO_POLICY" \
  --reward_threshold "$REWARD_THRESHOLD" \
  --output_file_path "$OUTPUT_FILE_PATH" \
  --results_path "$RESULTS_PATH" \
  --DATA_DIR "$DATA_DIR" \
  --variance "$VARIANCE" \
  --variance_min "$VARIANCE_MIN" \
  --variance_max "$VARIANCE_MAX" \
  --percentage "$PERCENTAGE" \
  --percentage_min "$PERCENTAGE_MIN" \
  --percentage_max "$PERCENTAGE_MAX" \
  --fxy_action "$FXY_ACTION" \
  --fxy_action_min "$FXY_ACTION_MIN" \
  --fxy_action_max "$FXY_ACTION_MAX" \
  --cxy_action "$CXY_ACTION" \
  --cxy_action_min "$CXY_ACTION_MIN" \
  --cxy_action_max "$CXY_ACTION_MAX"

  #!/bin/bash

TIMESTEPS=1000000
EPISODES=10
RESET_NUM_TIMESTEPS=False
ENT_COEF=0.01
LEARNING_RATE=0.0003
N_STEPS=2048
BATCH_SIZE=64
N_EPOCHS=10
GAMMA=0.99
CLIP_RANGE=0.20
MAX_GRAD_NORM=0.5
VERBOSE=1
RANDOM_SEED=42
DEVICE="cuda"
PPO_POLICY="MlpPolicy"
REWARD_THRESHOLD=0.05
OUTPUT_FILE_PATH="Results/state_reward_6-3_5.31.txt"
RESULTS_PATH="Results/results_6-3_5.31.csv"
DATA_DIR="/mnt/ST8000/jialei/SFM_E/data/6-3"
VARIANCE=0.1
VARIANCE_MIN=0.3
VARIANCE_MAX=2.1
PERCENTAGE=0.05
PERCENTAGE_MIN=0.0
PERCENTAGE_MAX=1.0
FXY_ACTION=0.02
FXY_ACTION_MIN=0.0
FXY_ACTION_MAX=1.0
CXY_ACTION=0.05
CXY_ACTION_MIN=0.0
CXY_ACTION_MAX=1.0

python main.py \
  --TIMESTEPS "$TIMESTEPS" \
  --episodes "$EPISODES" \
  --reset_num_timesteps "$RESET_NUM_TIMESTEPS" \
  --ent_coef "$ENT_COEF" \
  --learning_rate "$LEARNING_RATE" \
  --n_steps "$N_STEPS" \
  --batch_size "$BATCH_SIZE" \
  --n_epochs "$N_EPOCHS" \
  --gamma "$GAMMA" \
  --clip_range "$CLIP_RANGE" \
  --max_grad_norm "$MAX_GRAD_NORM" \
  --verbose "$VERBOSE" \
  --random_seed "$RANDOM_SEED" \
  --device "$DEVICE" \
  --ppo_Policy "$PPO_POLICY" \
  --reward_threshold "$REWARD_THRESHOLD" \
  --output_file_path "$OUTPUT_FILE_PATH" \
  --results_path "$RESULTS_PATH" \
  --DATA_DIR "$DATA_DIR" \
  --variance "$VARIANCE" \
  --variance_min "$VARIANCE_MIN" \
  --variance_max "$VARIANCE_MAX" \
  --percentage "$PERCENTAGE" \
  --percentage_min "$PERCENTAGE_MIN" \
  --percentage_max "$PERCENTAGE_MAX" \
  --fxy_action "$FXY_ACTION" \
  --fxy_action_min "$FXY_ACTION_MIN" \
  --fxy_action_max "$FXY_ACTION_MAX" \
  --cxy_action "$CXY_ACTION" \
  --cxy_action_min "$CXY_ACTION_MIN" \
  --cxy_action_max "$CXY_ACTION_MAX"