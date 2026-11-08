def collatz_function(n: int):
    print(f"Collatz Sequence for n={n}:")
    print(n)
    counts = []
    counts.append(n)
    while n != 1:
        if n % 2 == 0:
            iter = n/2
        else:
            iter = 3*n+1
        n = iter
        counts.append(n)
    return counts
# hold n sequences in a dictionary 
collatz_dict = {}
for i in range(100):
    collatz_dict[i+1] = collatz_function(i+1)

print(collatz_dict)