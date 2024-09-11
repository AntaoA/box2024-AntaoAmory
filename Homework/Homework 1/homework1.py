import random
import matplotlib.pyplot as plt


with open('Homework/Homework 1/genome_hw1.fa', 'r') as file:
    f = file.read()

list_seq = []
for line in f.split('>')[1:]:
    line = line.replace('\n', '')
    n = line.find("sequence")
    line = line[n+8:].strip()
    list_seq.append(line)

n_seq = random.randint(0, len(list_seq)-1)
random_seq = list_seq[n_seq]

print("Random sequence length: ", len(random_seq))

DNA = ['A', 'C', 'G', 'T']
random_DNA = "".join(random.choices(DNA, k = len(random_seq)))

def fibonacci_word(n):
    if n == 0:
        return ""
    if n == 1:
        return "A"
    sequence_start = "A"
    for _ in range(n-1):
        sequence = ""
        for i in range(len(sequence_start)):
            if sequence_start[i] == 'A':
                sequence = sequence + "AB"
            else:
                sequence = sequence + "A" 
        sequence_start = sequence 
    return sequence

def fibonacci_length(k):
    n = 0
    while len(fibonacci_word(n)) < k:
        n += 1
    return fibonacci_word(n)[:k]


fib = fibonacci_length(len(random_seq))

def F(seq, k):
    #compute the number of distinct subwords of length k in the sequence
    liste_subseq = [seq[i:i+k] for i in range(len(seq)-k+1)]
    ens_subseq = set(liste_subseq)
    return len(ens_subseq)


# Générer des valeurs entières pour k
k_values = range(1, 20)  # Plage de k de -10 à 10 (inclus)

# Calculer les valeurs de F(k)
F_values = [F(random_seq, k) for k in k_values]
F_values_DNA = [F(random_DNA, k) for k in k_values]
F_values_fib = [F(fibonacci_word(10), k) for k in k_values]
# Tracer le graphique
plt.figure(figsize=(10, 6))
plt.plot(k_values, F_values, 'o-', label='F(k)', color='blue')  # 'o-' pour les points et lignes
plt.plot(k_values, F_values_DNA, 'o-', label='F(k)', color='red')  # 'o-' pour les points et lignes
plt.plot(k_values, F_values_fib, 'o-', label='F(k)', color='green')  # 'o-' pour les points et lignes

plt.xlabel('k')
plt.ylabel('F(k)')
plt.title('Graphique de F(k) en fonction de k')
plt.legend()
plt.grid(True)
plt.show()