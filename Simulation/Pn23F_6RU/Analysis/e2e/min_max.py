# Open the file and read the contents
with open('Pn23F_6RU_V2_0_to_1000ns_e2e.txt', 'r') as file:
    # Read lines and extract the right column as floats
    values = [float(line.split()[1]) for line in file]

# Find the minimum and maximum values
min_value = min(values)
max_value = max(values)

# Print the results
print(f"Minimum value: {min_value}")
print(f"Maximum value: {max_value}")
