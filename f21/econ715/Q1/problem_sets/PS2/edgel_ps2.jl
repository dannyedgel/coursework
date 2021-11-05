#==
    This file completes all of the tasks for 1(e) of problem
    set 2 for Econ 715

    Date created:  04 November 2021
    Last modified: 04 November 2021
    Author: Danny Edgel
==#


using Distributed, SharedArrays

#add processes for parallelization
workers()
addprocs(4)