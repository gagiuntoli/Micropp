/*
 * This source code is part of MicroPP: a finite element library
 * to solve microstructural problems for composite materials.
 *
 * Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#ifndef INSTRUMENT_HPP
#define INSTRUMENT_HPP

#ifdef NDEBUG // For optimized compilation

#define INST_INIT
#define INST_START
#define INST_END
#define INST_FINAL

#else // For debug compilation

#include <unordered_map>
#include <vector>
#include <stack>
#include <string>
#include <atomic>

#include "util.hpp"

#define INST_INIT instrument::init()
#define INST_START instrument::registerStart(__FUNCTION__)
#define INST_END instrument::registerEnd(__FUNCTION__)
#define INST_FINAL instrument::finalize()

using namespace std;

typedef vector<struct timespec> timevect;
typedef stack<struct timespec> timestack;

class instrument {
	private:
		static unordered_map<string, timevect> times;                     // Final register time
		static thread_local  unordered_map<string, timestack> startTime;  // Calls stack
		static struct timespec initialTime;                         // Time collection
		static atomic<int> instances;

	public:
		static void init();

		static inline void registerStart(const string funct)
		{
			struct timespec time;
			getTime(time);
			startTime[funct].push(time);
		}

		static inline void registerEnd(const string funct)
		{
			struct timespec time2;
			getTime(time2);

			const struct timespec time1 = startTime[funct].top();
			startTime[funct].pop();

			times[funct].push_back(getDt(time1, time2));
		}

		static void finalize();
};

#endif // NDEBUG

#endif //INSTRUMENT_HPP

