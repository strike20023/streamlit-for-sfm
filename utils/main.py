from Env import SFMEnv
from stable_baselines3 import PPO
from stable_baselines3.common.logger import configure
import argparse
import sys
import os
import datetime



sys.path.append('./')

_PARSER = argparse.ArgumentParser('ppo_model')
_PARSER.add_argument('--TIMESTEPS', type=int, default=1000000, help="Deployment of each model learning")
_PARSER.add_argument('--episodes', type=int, default=40, help="Number of times the model was saved")
_PARSER.add_argument('--reset_num_timesteps', type=bool, default=False, help="Whether to start training from the current position")
_PARSER.add_argument('--ppo_run', type=str, default="ppo_run", help="Save the name of the log")
_PARSER.add_argument('--ent_coef', type=float, default=0.01, help="Entropy coefficient")
_PARSER.add_argument('--learning_rate', type=float, default=0.0003, help="Learning rate")
_PARSER.add_argument('--n_steps', type=int, default=2048, help="Number of steps")
_PARSER.add_argument('--batch_size', type=int, default=64, help="Batch size")
_PARSER.add_argument('--n_epochs', type=int, default=10, help="Number of epochs")
_PARSER.add_argument('--gamma', type=float, default=0.99, help="Discount factor")
_PARSER.add_argument('--clip_range', type=float, default=0.20, help="Clip range")
_PARSER.add_argument('--max_grad_norm', type=float, default=0.5, help="Maximum gradient")
_PARSER.add_argument('--use_sde', type=bool, default=False, help="Whether to use sde")
_PARSER.add_argument('--verbose', type=int, default=1, help="ppo verbose")
_PARSER.add_argument('--random_seed', type=int, default=42, help="random seed")
_PARSER.add_argument('--device', type=str, default="cuda", help="device")
_PARSER.add_argument('--ppo_Policy', type=str, default="MlpPolicy", help="ppo verbose")
_PARSER.add_argument('--reward_threshold', type=float, default=0.05, help="reward_threshold")
_PARSER.add_argument('--output_file_path', type=str, default="Results/state_reward.txt", help="output_file_path")
_PARSER.add_argument('--results_path', type=str, default="Results/results.csv", help="results_path")
_PARSER.add_argument('--DATA_DIR', type=str, default="/mnt/ST8000/zhenhanbai/RL/SFM", help="image dir")
_PARSER.add_argument('--variance', type=float, default=0.1, help="Increase the step size of the variance")
_PARSER.add_argument('--variance_min', type=float, default=0.3, help="the min scale of variance")
_PARSER.add_argument('--variance_max', type=float, default=2.1, help="the max scale of variance")
_PARSER.add_argument('--percentage', type=float, default=0.05, help="Increase the step size of the percentage")
_PARSER.add_argument('--percentage_min', type=float, default=0.0, help="the min scale of percentage")
_PARSER.add_argument('--percentage_max', type=float, default=1.0, help="the max scale of percentage")
_PARSER.add_argument('--fxy_action', type=float, default=0.01, help="Increase the step size of the focal length error")
_PARSER.add_argument('--fxy_action_min', type=float, default=0.0, help="the min scale of focal length")
_PARSER.add_argument('--fxy_action_max', type=float, default=1.0, help="the max scale of focal length")
_PARSER.add_argument('--cxy_action', type=float, default=0.05, help="Increase the step size of the Light shaft error")
_PARSER.add_argument('--cxy_action_min', type=float, default=0.0, help="the min scale of Light shaft")
_PARSER.add_argument('--cxy_action_max', type=float, default=1.0, help="the max scale of Light shaft")

_ARGS = _PARSER.parse_args()

print(_ARGS)

dirs_to_create = ["Results/Training_models", "Results/Logs", "Results"]
[os.makedirs(dir) for dir in dirs_to_create if not os.path.exists(dir)]
new_logger = configure(dirs_to_create[1], ["stdout", "csv", "tensorboard"])

if __name__ == "__main__":
    start_time = datetime.datetime.now()
    env = SFMEnv(_ARGS)
    env.reset()
    model = PPO('MlpPolicy', 
                env, 
                ent_coef=_ARGS.ent_coef,
                learning_rate=_ARGS.learning_rate,
                n_steps=_ARGS.n_steps,
                batch_size=_ARGS.batch_size,
                n_epochs=_ARGS.n_epochs,
                gamma=_ARGS.gamma,
                clip_range=_ARGS.clip_range,
                max_grad_norm=_ARGS.max_grad_norm,
                use_sde=_ARGS.use_sde,
                verbose=_ARGS.verbose,
                seed=_ARGS.random_seed,
                device=_ARGS.device)
    model.set_logger(new_logger)

    for i in range(_ARGS.episodes):
        model.save(f"{dirs_to_create[0]} / {_ARGS.TIMESTEPS * i}")
        model.learn(total_timesteps=_ARGS.TIMESTEPS, reset_num_timesteps=_ARGS.reset_num_timesteps, tb_log_name=_ARGS.ppo_run)
        # 尝试 episodes 写入
        env.Calculate_state_results()

    end_time = datetime.datetime.now()
    print("----------start_time----------\n", start_time)
    print("----------end_time----------\n", end_time)
    print("--------------------\n", end_time - start_time)