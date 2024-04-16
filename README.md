# Calculul paralel al unei expresii

Folosind biblioteca OpenMPI, am implementat un algoritm de calcul paralelizat al unei fractii, rezultatul fiind salvat in fisier.

Numarul minim de procese pe care se poate rula programul este 4, numarul dat de procese trebuie neaparat sa fie par.

Datele sunt distribuite catre procesele worker de catre procesele MASTER din fiecare comunicator (lowComm si highComm). Procesele master NU participa si ele la calcul, ci doar aduna rezultatele de la workeri si le trimite mai departe catre procesul MASTER al programului pentru calculul final si scrierea in fisierul de iesire.
