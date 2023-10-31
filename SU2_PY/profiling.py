#!/usr/bin/env python

## \file profiling.py
#  \brief Python script for postprocessing the SU2 custom profiling (profiling.csv)
#  \author T. Economon
#  \version 8.0.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

from optparse import OptionParser
from pylab import *
from numpy import *
from matplotlib import pyplot as plt
from matplotlib import mlab

parser = OptionParser()
parser.add_option(
    "-f", "--file", dest="file", help="profiling CSV file", metavar="FILE"
)
(options, args) = parser.parse_args()

# Store the file name
filename = options.file

# Load the csv file with the profiling data
profile = mlab.csv2rec(filename, comments="#", skiprows=0, checkrows=0)

# Get total number of groups
maxID = 0
for val in range(len(profile.function_name)):
    if profile.function_id[val] > maxID:
        maxID = profile.function_id[val]

# Get some arrays for sorting out the groups
labels = [[] for i in range(maxID)]
fracs = [[] for i in range(maxID)]
explode = [[] for i in range(maxID)]
calls = [[] for i in range(maxID)]

# Process the profiling data into group IDs
for val in range(len(profile.function_name)):
    labels[profile.function_id[val] - 1].append(profile.function_name[val])
    fracs[profile.function_id[val] - 1].append(profile.avg_total_time[val])
    explode[profile.function_id[val] - 1].append(0)
    calls[profile.function_id[val] - 1].append(profile.n_calls[val])

# Loop over each of the group IDs and make figures
for val in range(maxID):

    # Create a Pie chart to see the time spent in each subroutine
    fig = plt.figure(figsize=[18, 8])
    ax = fig.add_subplot(121)
    ax.set_title("Total Time Spent in Each Function")

    # Sort the pieces to make it pretty
    fracs[val], labels[val], calls[val] = (
        list(x) for x in zip(*sorted(zip(fracs[val], labels[val], calls[val])))
    )

    # Call to make the pie chart
    pie_wedge_collection = ax.pie(
        fracs[val],
        explode=explode[val],
        labels=labels[val],
        labeldistance=1.05,
        autopct="%1.1f%%",
        shadow=False,
        startangle=0,
    )
    for pie_wedge in pie_wedge_collection[0]:
        pie_wedge.set_edgecolor("white")

    # Sort the number of calls again for the bar chart
    calls[val], labels[val] = (
        list(x) for x in zip(*sorted(zip(calls[val], labels[val])))
    )

    # Create a bar chart for the number of function calls
    ax = fig.add_subplot(122)
    ax.set_title("Number of Function Calls")
    width = 0.35
    ax.bar(range(len(calls[val])), calls[val], width=width)
    ax.set_xticks(np.arange(len(calls[val])) + width / 2)
    ax.set_xticklabels(labels[val])
    ax.set_xlabel("Function")
    ax.set_ylabel("Calls")
    fig.autofmt_xdate()
    fig.subplots_adjust(wspace=0.5)

    # Save a figure for this group
    filename = "profile_group_" + str(val) + ".png"
    fig.savefig(filename, format="png")

    # Uncomment the next line to open the plots on the screen
    show()
