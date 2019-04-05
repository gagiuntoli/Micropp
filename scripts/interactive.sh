#!/bin/bash
# @ job_name = interactive
# @ initialdir = ./
# @ output = ./interactive_%j.out
# @ error = ./interactive_%j.err
# @ wall_clock_limit = 01:00:00
# @ tasks_per_node = 1
# @ cpus_per_task = 16
# @ gpus_per_node = 2
# @ features = k80
# @ X11 = 1
sleep 2h
