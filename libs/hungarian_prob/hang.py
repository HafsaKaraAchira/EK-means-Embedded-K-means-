# Python program for the above approach
import dlib

# Function to find out the best
# assignment of people to jobs so that
# total cost of the assignment is minimized
def minCost(arr):

	# Call the max_cost_assignment() function
	# and store the assignment
	assignment = dlib.max_cost_assignment(arr)

	# Print the optimal cost
	print(dlib.assignment_cost(arr, assignment))

# Driver Code

# Given 2D array
arr = dlib.matrix([[3, 5], [10, 1]])

# Function Call
minCost(arr)
