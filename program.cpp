#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>

#define N 200
#define TAG_VECTOR_X 1
#define TAG_VECTOR_Y 2
#define TAG_MATRIX_A 3
#define TAG_SUMA_JOS 4
#define TAG_SUMA_SUS 5
#define inDebugMode 0
#define showReceivedData 0
#define MASTER_PROC_ID 0

using namespace std;
// serializam o bucata din matrice in vector unidimensional, pornind de la un anumit rand (dat de offset) si cate randuri trimitem
void serialize(const vector<vector<float>> &data, float *buffer, int offset, int rowsToSend)
{
    int index = 0;
    for (int i = offset; i < offset + rowsToSend; ++i)
    {
        for (float value : data[i])
        {
            buffer[index++] = value;
        }
    }
}

// deserializam o bucata din matrice din vector unidimensional in vector bidimensional, cunoscand numarul de randuri primite
void deserialize(const float *buffer, vector<vector<float>> &data, int rowsToReceive)
{
    data.clear();
    int index = 0;
    for (int i = 0; i < rowsToReceive; ++i)
    {
        vector<float> row;
        for (int j = 0; j < N; ++j)
        {
            row.push_back(buffer[index++]);
        }
        data.push_back(row);
    }
}

float calculateLow(vector<float> xData, vector<float> yData, int numOfElements)
{
    float suma = 0.0f;
    for (int i = 0; i < numOfElements; i++)
    {
        suma += xData[i] * yData[i];
    }
    return suma;
}

float calculateHigh(vector<float> xData, vector<float> yData, vector<vector<float>> matData, int noOfElements)
{
    float suma = 0.0f;
    for (int i = 0; i < noOfElements; i++)
    {
        for (int j = 0; j < N; j++)
        {
            suma += xData[i] * matData[i][j] * yData[j];
        }
    }
    return suma;
}

int main(int argc, char **argv)
{
    int mainRank, nrProc, secondRank;
    MPI_Comm lowComm, highComm;

    // init mpi
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mainRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nrProc);

    if (nrProc % 2 != 0 || nrProc < 4)
    {
        cout << "Numarul de procese trebuie sa fie par si mai mare ca 2!" << endl;
        MPI_Finalize();
        return 0;
    } // verificam daca nr de procese este par inainte de a face alocari
    // verificam si ca numarul minim de procese sa fie 4 deoarece neavand voie ca procesul master pe fiecare subgrup
    // sa faca calcule ci doar sa adune si sa trimita date mai departe, am ajunge la o eroare unde nu avem procese worker.

    int lowGroupMaxRank = nrProc / 2;
    int highGroupMinRank = nrProc - lowGroupMaxRank;

    if (mainRank < lowGroupMaxRank)
    {
        MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &lowComm);
        MPI_Comm_rank(lowComm, &secondRank);
    }
    else
    {
        MPI_Comm_split(MPI_COMM_WORLD, 1, 0, &highComm);
        MPI_Comm_rank(highComm, &secondRank);
    }

    if (mainRank < lowGroupMaxRank)
    {
        if (inDebugMode == 1)
        {
            cout << "Proces din grupul LOW cu rank principal: " << mainRank << " si rank in grupul LOW: " << secondRank << endl;
        }

        int masterLow = 0;
        int elementsPerProcess = N / (lowGroupMaxRank - 1); // calculam numarul de elemente care va fi procesat de fiecare worker
        int r = N % (lowGroupMaxRank - 1);
        float localLowSum = 0.0f;
        float totalLowSum = 0.0f;
        if (mainRank == lowGroupMaxRank - 1)
        {
            elementsPerProcess += r;
        } // daca nu se poate imparti egal, ultimul proces va primi restul de elemente

        if (inDebugMode == 1)
        {
            if (secondRank != 0)
            {
                cout << "Procesul cu ID: " << mainRank << " va primi " << elementsPerProcess << " date."
                     << " WORKER LOW GROUP" << endl;
            }
        }

        if (secondRank == 0 && mainRank < lowGroupMaxRank)
        {
            int masterLowGroup = mainRank;
            vector<float> xData;
            vector<float> yData;

            ifstream xFile("x.dat");
            if (!xFile.is_open())
            {
                cerr << "Eroare deschidere x.dat." << endl;
                MPI_Finalize();
                return 0;
            }
            float value;
            while (xFile >> value)
            {
                xData.push_back(value);
            }
            xFile.close();

            ifstream yFile("y.dat");
            if (!yFile.is_open())
            {
                cerr << "Eroare deschidere y.dat." << endl;
                MPI_Finalize();
                return 0;
            }
            while (yFile >> value)
            {
                yData.push_back(value);
            }
            yFile.close();

            // trimitem date la fiecare worker
            for (int destRank = 1; destRank < lowGroupMaxRank; ++destRank) // pornim de la 1 ca sa nu trimitem si procesului master
            {
                // calculam offset-ul de date pentru procesul curent la care trimitem
                int offset = (destRank - 1) * elementsPerProcess;

                // calculam numarul de elemente care trebuie trimise
                int elementsToSend = elementsPerProcess;
                if (destRank == lowGroupMaxRank - 1) // daca numarul de workeri este impar ultimul proces primeste si restul de date
                {
                    elementsToSend += r;
                }

                MPI_Send(xData.data() + offset, elementsToSend, MPI_FLOAT, destRank, TAG_VECTOR_X, lowComm);
                MPI_Send(yData.data() + offset, elementsToSend, MPI_FLOAT, destRank, TAG_VECTOR_Y, lowComm);
            }
        }
        else
        {
            vector<float> localXData(elementsPerProcess);
            vector<float> localYData(elementsPerProcess);

            MPI_Recv(localXData.data(), elementsPerProcess, MPI_FLOAT, masterLow, TAG_VECTOR_X, lowComm, MPI_STATUS_IGNORE);
            MPI_Recv(localYData.data(), elementsPerProcess, MPI_FLOAT, masterLow, TAG_VECTOR_Y, lowComm, MPI_STATUS_IGNORE);

            if (inDebugMode == 1 && showReceivedData == 1)
            {
                for (int i = 0; i < elementsPerProcess; ++i)
                {
                    cout << "xData[" << i << "]: " << localXData[i] << ", yData[" << i << "]: " << localYData[i] << endl;
                }
            }

            localLowSum = calculateLow(localXData, localYData, elementsPerProcess);
        }
        MPI_Reduce(&localLowSum, &totalLowSum, 1, MPI_FLOAT, MPI_SUM, masterLow, lowComm);
        if (secondRank == 0)
        {
            if (inDebugMode == 1)
            {
                cout << "Suma totala grup LOW: " << totalLowSum << endl;
            }

            MPI_Send(&totalLowSum, 1, MPI_FLOAT, MASTER_PROC_ID, TAG_SUMA_JOS, MPI_COMM_WORLD);
        }
    }
    else
    {
        if (mainRank >= highGroupMinRank)
        {
            if (inDebugMode == 1)
            {
                cout << "Proces din grupul HIGH cu rank principal: " << mainRank << " si rank din grupul HIGH: " << secondRank << endl;
            }

            int masterHigh = 0;
            int elementsPerProcess = N / ((nrProc - highGroupMinRank) - 1); // calculam numarul de elemente care va fi procesat de fiecare worker
            int r = N % ((nrProc - highGroupMinRank) - 1);
            float localSumHigh = 0.0f;
            float totalSumHigh = 0.0f;
            if (mainRank == (nrProc - 1))
            {
                elementsPerProcess += r;
            } // daca nu se poate imparti egal, ultimul proces va primi restul de elemente

            if (inDebugMode == 1)
            {
                if (secondRank != 0)
                {
                    cout << "Procesul cu ID: " << mainRank << " va primi " << elementsPerProcess << " date."
                         << " WORKER HIGH GROUP" << endl;
                }
            }
            if (secondRank == 0 && mainRank >= highGroupMinRank)
            {
                vector<float> xData;
                vector<float> yData;
                vector<vector<float>> matData(N, vector<float>(N));

                ifstream xFile("x.dat");
                if (!xFile.is_open())
                {
                    cerr << "Eroare deschidere x.dat." << endl;
                    MPI_Finalize();
                    return 0;
                }
                float value;
                while (xFile >> value)
                {
                    xData.push_back(value);
                }
                xFile.close();

                ifstream yFile("y.dat");
                if (!yFile.is_open())
                {
                    cerr << "Eroare deschidere y.dat." << endl;
                    MPI_Finalize();
                    return 0;
                }
                while (yFile >> value)
                {
                    yData.push_back(value);
                }
                yFile.close();

                ifstream matFile("mat.dat");
                if (!matFile.is_open())
                {
                    cerr << "Eroare deschidere mat.dat." << endl;
                    MPI_Finalize();
                    return 0;
                }

                for (int i = 0; i < N; ++i)
                {
                    for (int j = 0; j < N; ++j)
                    {
                        matFile >> matData[i][j];
                    }
                }
                matFile.close();

                for (int destRank = 1; destRank < (nrProc - highGroupMinRank); ++destRank)
                {
                    int offset = (destRank - 1) * elementsPerProcess;
                    int elementsToSend = elementsPerProcess;
                    if (destRank == (nrProc - highGroupMinRank) - 1) //la fel ca la cele din LOW
                    {
                        elementsToSend += r;
                    }

                    MPI_Send(xData.data() + offset, elementsToSend, MPI_FLOAT, destRank, TAG_VECTOR_X, highComm);
                    MPI_Send(yData.data(), N, MPI_FLOAT, destRank, TAG_VECTOR_Y, highComm);

                    // serializam datele matricei corespunzator cu offset si numarul de linii de trimis
                    auto *buffer = new float[elementsToSend * N];
                    serialize(matData, buffer, offset, elementsToSend);
                    MPI_Send(buffer, elementsToSend * N, MPI_FLOAT, destRank, TAG_MATRIX_A, highComm);
                    delete[] buffer;
                }
            }
            else
            {
                vector<float> localXData(elementsPerProcess);
                vector<float> localYData(N);
                vector<vector<float>> localMatData(elementsPerProcess, vector<float>(N));

                MPI_Recv(localXData.data(), elementsPerProcess, MPI_FLOAT, masterHigh, TAG_VECTOR_X, highComm, MPI_STATUS_IGNORE);
                MPI_Recv(localYData.data(), N, MPI_FLOAT, masterHigh, TAG_VECTOR_Y, highComm, MPI_STATUS_IGNORE);

                auto *buffer = new float[elementsPerProcess * N];
                MPI_Recv(buffer, elementsPerProcess * N, MPI_FLOAT, masterHigh, TAG_MATRIX_A, highComm, MPI_STATUS_IGNORE);

                // deserializam matricea corespunzator
                deserialize(buffer, localMatData, elementsPerProcess);
                delete[] buffer;
                if (inDebugMode == 1 && showReceivedData == 1)
                {
                    cout << "Procesul " << mainRank << " a primit datele din vectorul X:" << endl;
                    for (int i = 0; i < elementsPerProcess; ++i)
                    {
                        cout << localXData[i] << " ";
                    }
                    cout << endl;

                    cout << "Procesul " << mainRank << " a primit datele din vectorul Y:" << endl;
                    for (int i = 0; i < N; ++i)
                    {
                        cout << localYData[i] << " ";
                    }
                    cout << endl;

                    cout << "Procesul " << mainRank << " a primit datele din matrice:" << endl;
                    for (int i = 0; i < elementsPerProcess; ++i)
                    {
                        for (int j = 0; j < N; ++j)
                        {
                            cout << localMatData[i][j] << " ";
                        }
                        cout << endl;
                    }
                }

                localSumHigh = calculateHigh(localXData, localYData, localMatData, elementsPerProcess);
            }
            MPI_Reduce(&localSumHigh, &totalSumHigh, 1, MPI_FLOAT, MPI_SUM, masterHigh, highComm);
            if (secondRank == 0)
            {
                if (inDebugMode == 1)
                {
                    cout << "Suma totala grup HIGH: " << totalSumHigh << endl;
                }
                MPI_Send(&totalSumHigh, 1, MPI_FLOAT, MASTER_PROC_ID, TAG_SUMA_SUS, MPI_COMM_WORLD);
            }
        }
    }
    if (mainRank == 0)
    {
        float rezultatJos = 0.0f, rezultatSus = 0.0f, cantitateTotala = 0.0f;
        MPI_Recv(&rezultatJos, 1, MPI_FLOAT, MPI_ANY_SOURCE, TAG_SUMA_JOS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rezultatSus, 1, MPI_FLOAT, MPI_ANY_SOURCE, TAG_SUMA_SUS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cantitateTotala = rezultatSus / rezultatJos;
        if (inDebugMode == 1)
        {
            cout << "AVG: " << cantitateTotala << endl;
        }

        ofstream outputFile("avg.txt");
        if (outputFile.is_open())
        {
            outputFile << "AVG= " << cantitateTotala << endl;
            outputFile.close();
            cout << "Rezultat salvat in fisierul avg.txt" << endl;
        }
        else
        {
            cerr << "Eroare la scrierea lui avg.txt." << endl;
        }
    }
    MPI_Finalize();

    return 0;
}
