import random
import time

n = 10000
list = [random.randint(1, 100) for x in range(n)]

print(list)

def merge(a, b):

    # Merges 2 Sorted Lists

    # i, j are counters for lists a, b
    list = []
    i = 0
    j = 0

    # Continue until counter reaches end of one of the input lists
    while i < len(a) and j < len(b):
        if a[i] <= b[j]:
            list.append(a[i])
            i += 1
        else:
            list.append(b[j])
            j += 1
    # Append end of remaining list
    if i < len(a):
        [list.append(_) for _ in a[i:]]
    else:
        [list.append(_) for _ in b[j:]]

    return list

def merge_sort(a):

    n = len(a)
    if n > 1:
        a = merge(merge_sort(a[:n//2]), merge_sort(a[n//2:]))

    return a

def selection_sort(a):

    for i in range(len(a)):
        min_idx = i
        for j in range(i + 1, len(a)):
            if a[min_idx] > a[j]:
                min_idx = j
        a[i], a[min_idx] = a[min_idx], a[i]

    return a

time1 = time.time()
print(merge_sort(list))
print(f'Merge sort on {n:,} size list: {time.time()-time1:.2f} seconds')

time1 = time.time()
print(selection_sort(list))
print(f'Selection sort on {n:,} size list: {time.time()-time1:.2f} seconds')




