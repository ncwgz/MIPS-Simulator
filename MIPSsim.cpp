/*
 * On my honor, I have neither given nor received unauthorized aid on this assignment
 *
 * Author:  Guozhi Wang
 * UFID:    9357-8947
 * Email:   wangguozhi@ufl.edu
 * Date:    2020-11-19
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <iomanip>

using namespace std;

struct Instruction {
    string assembly;
    int category;
    int opcode;
    string operation;
    unsigned int r1;
    unsigned int r2;
    unsigned int r3;
    long lValue;
    int sValue;
    int index;
};

class Processor {
private:
    unsigned int PC;
    long* registers;
    map<int, long>* memory;
    map<int, Instruction>* instructions;
    vector<int> *alteredMemory;
    int cycleCounter;

public:
    Processor () {
        this->PC = 256;
        this->registers = new long[32];
        this->memory = new map<int, long>();
        this->instructions = new map<int, Instruction>();
        this->alteredMemory = new vector<int>();

        for (int i = 0; i < 32; ++i) {
            registers[i] = 0;
        }

        this->cycleCounter = 1;
    }

    void addInstruction(Instruction i) {
        this->instructions->insert(pair<int, Instruction>(this->PC, i));
        PC += 4;
    }

    string getAssembly(Instruction i) {
        stringstream ss;
        ss << i.index << "\t";
        if (i.category == 0) {
            ss << i.sValue << endl;
        }
        else if (i.category == 1) {
            ss << i.operation << " ";
            switch (i.opcode) {
                case 0: { // J
                    ss << "#" << i.lValue << endl;
                    break;
                }
                case 1: { // JR
                    ss << "R" << i.r1 << endl;
                    break;
                }
                case 2: { // BEQ
                    ss << "R" << i.r1 << ", ";
                    ss << "R" << i.r2 << ", ";
                    ss << "#" << i.lValue << endl;
                    break;
                }
                case 3: { // BLGZ
                    // the same to case 4(BGTZ)
                }
                case 4: { // BGTZ
                    ss << "R" << i.r1 << ", ";
                    ss << "#" << i.lValue << endl;
                    break;
                }
                case 5: { // BREAK
                    ss << endl;
                    break;
                }
                case 6: { // SW
                    // the same to case 7(LW)
                }
                case 7: { // LW
                    ss << "R" << i.r2 << ", ";
                    ss << "" << i.lValue;
                    ss << "(R" << i.r1 << ")" << endl;
                    break;
                }
                case 8: { // SLL
                    // the same to case 10(SRA)
                }
                case 9: { // SRL
                    // the same to case 10(SRA)
                }
                case 10: { // SRA
                    ss << "R" << i.r2 << ", ";
                    ss << "R" << i.r1 << ", ";
                    ss << "#" << i.lValue << endl;
                    break;
                }
                case 11: { // NOP
                    ss << endl;
                    break;
                };
            }
        }
        else if (i.category == 3) {
            ss << i.operation << " ";
            switch (i.opcode) {
                case 0: { // ADD
                    // the same to case 7(SLT)
                }
                case 1: { // SUB
                    // the same to case 7(SLT)
                }
                case 2: { // MUL
                    // the same to case 7(SLT)
                }
                case 3: { // AND
                    // the same to case 7(SLT)
                }
                case 4: { // OR
                    // the same to case 7(SLT)
                }
                case 5: { // XOR
                    // the same to case 7(SLT)
                }
                case 6: { // NOR
                    // the same to case 7(SLT)
                }
                case 7: { // SLT
                    ss << "R" << i.r3 << ", ";
                    ss << "R" << i.r1 << ", ";
                    ss << "R" << i.r2 << endl;
                    break;
                }
                case 8: { // ADDI
                    ss << "R" << i.r2 << ", ";
                    ss << "R" << i.r1 << ", ";
                    ss << "#" << i.sValue << endl;
                    break;
                }
                case 9: { // ANDI
                    // the same to case 11(XORI)
                }
                case 10: { // ORI
                    // the same to case 11(XORI)
                }
                case 11: { // XORI
                    ss << "R" << i.r2 << ", ";
                    ss << "R" << i.r1 << ", ";
                    ss << "#" << i.lValue << endl;
                    break;
                }
            }
        }
        return ss.str();
    }

    void writeDisassemblyFile() {
        ofstream file("disassembly.txt", ios::out);
        int tempCounter = 256;
        for (auto imap : *this->instructions) {
            Instruction i = imap.second;
            file << i.assembly << "\t";
            file << getAssembly(i);
            tempCounter += 4;
        }

        file.close();
    }

    void loadProgramData() {
        for (auto imap : *this->instructions) {
            Instruction i = imap.second;
            if (i.category == 0) {
                this->memory->insert(pair<int, long>(i.index, i.sValue));
            }
        }
    }

    void writeSimulationFile(ofstream& file, Instruction i) {
        file << "--------------------" << endl;
        file << "Cycle " << cycleCounter << ":\t";
        file << getAssembly(i) << endl;

        file << "Registers" << endl;
        for (int r = 0; r < 32; ++r) {
            if (r % 8 == 0) {
                file << "R";
                file << setfill('0') << setw(2) << r;
                file << ":";
            }
            file << "\t" << registers[r];
            if (r % 8 == 7) {
                file << endl;
            }
        }
        file << endl;

        file << "Data" << endl;
        int counter = 0;
        for (auto m : *memory) {
            if (counter == 0) {
                file << m.first << ":";
            }
            file << "\t" << m.second;
            counter ++;
            if (counter == 8) {
                file << endl;
                counter = 0;
            }
        }
    }

    void simulate() {
        this->cycleCounter = 0;
        this->alteredMemory->clear();
        bool isToBreak = false;
        this->PC = 256;
        ofstream file("simulation.txt", ios::out);
        while (true) {
            cycleCounter ++;
            Instruction i = instructions->find(PC)->second;
            PC += 4;
            if (i.category == 1) {
                switch (i.opcode) {
                    case 0: { // J
                        PC = i.lValue;
                        break;
                    }
                    case 1: { // JR
                        PC = registers[i.r1];
                        break;
                    }
                    case 2: { // BEQ
                        if (registers[i.r1] == registers[i.r2]) {
                            PC += i.lValue;
                        }
                        break;
                    }
                    case 3: { // BLTZ
                        if (registers[i.r1] < 0) {
                            PC += i.lValue;
                        }
                        break;
                    }
                    case 4: { // BGTZ
                        if (registers[i.r1] > 0) {
                            PC += i.lValue;
                        }
                        break;
                    }
                    case 5: { // BREAK
                        isToBreak = true;
                        break;
                    }
                    case 6: { // SW
                        memory->find(registers[i.r1] + i.lValue)->second = registers[i.r2];
                        break;
                    }
                    case 7: { // LW
                        registers[i.r2] = memory->find(registers[i.r1] + i.lValue)->second;
                        break;
                    }
                    case 8: { // SLL
                        registers[i.r2] = registers[i.r1] << i.lValue;
                        break;
                    }
                    case 9: { // SRL
                        registers[i.r2] = ((unsigned int)registers[i.r1]) >> i.lValue;
                        break;
                    }
                    case 10: { // SRA
                        registers[i.r2] = registers[i.r1] >> i.lValue;
                        break;
                    }
                    case 11: { // NOP
                        break;
                    }
                }
            }
            else if (i.category == 3) {
                switch (i.opcode) {
                    case 0: { // ADD
                        // rd = rs op rt
                        registers[i.r3] = registers[i.r1] + registers[i.r2];
                        break;
                    }
                    case 1: { // SUB
                        registers[i.r3] = registers[i.r1] - registers[i.r2];
                        break;
                    }
                    case 2: { // MUL
                        registers[i.r3] = registers[i.r1] * registers[i.r2];
                        break;
                    }
                    case 3: { // AND
                        registers[i.r3] = registers[i.r1] & registers[i.r2];
                        break;
                    }
                    case 4: { // OR
                        registers[i.r3] = registers[i.r1] | registers[i.r2];
                        break;
                    }
                    case 5: { // XOR
                        registers[i.r3] = registers[i.r1] ^ registers[i.r2];
                        break;
                    }
                    case 6: { // NOR
                        registers[i.r3] = ~(registers[i.r1] | registers[i.r2]);
                        break;
                    }
                    case 7: { // SLT
                        // rd = rs < rt
                        if (registers[i.r1] < registers[i.r2]) {
                            registers[i.r3] = 1;
                        } else {
                            registers[i.r3] = 0;
                        }
                        break;
                    }
                    case 8: { // ADDI
                        // rt = rs op immediate
                        registers[i.r2] = registers[i.r1] + i.sValue;
                        break;
                    }
                    case 9: { // ANDI
                        registers[i.r2] = registers[i.r1] & i.lValue;
                        break;
                    }
                    case 10: { // ORI
                        registers[i.r2] = registers[i.r1] | i.lValue;
                        break;
                    }
                    case 11: { // XORI
                        registers[i.r2] = registers[i.r1] ^ i.lValue;
                        break;
                    }
                }
            }

            writeSimulationFile(file, i);
            if (isToBreak) break;
        }
        file << endl;
        file.close();
    }
};

int main (int argc, char** argv) {

    if (argc != 2) {
        cout << "Please provide the file path." << endl;
        return 0;
    }

    Processor* processor = new Processor();
    int index = 256;

    string filepath = argv[1];

    string instruction;
    ifstream file(filepath, ios::in);
    bool isBreak = false;
    while (file >> instruction) {
        int category = stoi(instruction.substr(0, 2), nullptr, 2);
        int opcode = stoi(instruction.substr(2, 4), nullptr, 2);
        Instruction i;
        i.index = index;
        index += 4;
        if (isBreak) {
            i.category = 0;
            i.sValue = stol(instruction, nullptr, 2);
            i.assembly = instruction;
            processor->addInstruction(i);
        } else {
            if (category == 1) {
                i.category = 1;
                i.opcode = opcode;
                i.assembly = instruction;
                switch (opcode) {
                    case 0: {
                        i.operation = "J";
                        long instr_index = stol(instruction.substr(6, 26), nullptr, 2) << 2;
                        i.lValue = instr_index;
                        processor->addInstruction(i);
                        break;
                    }
                    case 1: {
                        i.operation = "JR";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        i.r1 = rs;
                        processor->addInstruction(i);
                        break;
                    }
                    case 2: {
                        i.operation = "BEQ";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        long offset = stol(instruction.substr(16, 16), nullptr, 2) << 2;
                        i.r1 = rs;
                        i.r2 = rt;
                        i.lValue = offset;
                        processor->addInstruction(i);
                        break;
                    }
                    case 3: {
                        i.operation = "BLTZ";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        long offset = stol(instruction.substr(16, 16), nullptr, 2) << 2;
                        i.r1 = rs;
                        i.lValue = offset;
                        processor->addInstruction(i);
                        break;
                    }
                    case 4: {
                        i.operation = "BGTZ";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        long offset = stol(instruction.substr(16, 16), nullptr, 2) << 2;
                        i.r1 = rs;
                        i.lValue = offset;
                        processor->addInstruction(i);
                        break;
                    }
                    case 5: {
                        i.operation = "BREAK";
                        processor->addInstruction(i);
                        isBreak = true;
                        break;
                    }
                    case 6: {
                        i.operation = "SW";
                        unsigned int base = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        long offset = stol(instruction.substr(16, 16), nullptr, 2);
                        i.r1 = base;
                        i.r2 = rt;
                        i.lValue = offset;
                        processor->addInstruction(i);
                        break;
                    }
                    case 7: {
                        i.operation = "LW";
                        unsigned int base = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        long offset = stol(instruction.substr(16, 16), nullptr, 2);
                        i.r1 = base;
                        i.r2 = rt;
                        i.lValue = offset;
                        processor->addInstruction(i);
                        break;
                    }
                    case 8: {
                        i.operation = "SLL";
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        unsigned int rd = stoul(instruction.substr(16, 5), nullptr, 2);
                        long sa = stol(instruction.substr(21, 5), nullptr, 2);
                        i.r1 = rt;
                        i.r2 = rd;
                        i.lValue = sa;
                        processor->addInstruction(i);
                        break;
                    }
                    case 9: {
                        i.operation = "SRL";
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        unsigned int rd = stoul(instruction.substr(16, 5), nullptr, 2);
                        long sa = stol(instruction.substr(21, 5), nullptr, 2);
                        i.r1 = rt;
                        i.r2 = rd;
                        i.lValue = sa;
                        processor->addInstruction(i);
                        break;
                    }
                    case 10: {
                        i.operation = "SRA";
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        unsigned int rd = stoul(instruction.substr(16, 5), nullptr, 2);
                        long sa = stol(instruction.substr(21, 5), nullptr, 2);
                        i.r1 = rt;
                        i.r2 = rd;
                        i.lValue = sa;
                        processor->addInstruction(i);
                        break;
                    }
                    case 11: {
                        i.operation = "NOP";
                        processor->addInstruction(i);
                        break;
                    }
                }
                continue;
            }
            else if (category == 3) {
                i.category = 3;
                i.opcode = opcode;
                i.assembly = instruction;
                switch (opcode) {
                    case 0: {
                        i.operation = "ADD";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        unsigned int rd = stoul(instruction.substr(16, 5), nullptr, 2);
                        i.r1 = rs;
                        i.r2 = rt;
                        i.r3 = rd;
                        processor->addInstruction(i);
                        break;
                    }
                    case 1: {
                        i.operation = "SUB";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        unsigned int rd = stoul(instruction.substr(16, 5), nullptr, 2);
                        i.r1 = rs;
                        i.r2 = rt;
                        i.r3 = rd;
                        processor->addInstruction(i);
                        break;
                    }
                    case 2: {
                        i.operation = "MUL";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        unsigned int rd = stoul(instruction.substr(16, 5), nullptr, 2);
                        i.r1 = rs;
                        i.r2 = rt;
                        i.r3 = rd;
                        processor->addInstruction(i);
                        break;
                    }
                    case 3: {
                        i.operation = "AND";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        unsigned int rd = stoul(instruction.substr(16, 5), nullptr, 2);
                        i.r1 = rs;
                        i.r2 = rt;
                        i.r3 = rd;
                        processor->addInstruction(i);
                        break;
                    }
                    case 4: {
                        i.operation = "OR";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        unsigned int rd = stoul(instruction.substr(16, 5), nullptr, 2);
                        i.r1 = rs;
                        i.r2 = rt;
                        i.r3 = rd;
                        processor->addInstruction(i);
                        break;
                    }
                    case 5: {
                        i.operation = "XOR";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        unsigned int rd = stoul(instruction.substr(16, 5), nullptr, 2);
                        i.r1 = rs;
                        i.r2 = rt;
                        i.r3 = rd;
                        processor->addInstruction(i);
                        break;
                    }
                    case 6: {
                        i.operation = "NOR";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        unsigned int rd = stoul(instruction.substr(16, 5), nullptr, 2);
                        i.r1 = rs;
                        i.r2 = rt;
                        i.r3 = rd;
                        processor->addInstruction(i);
                        break;
                    }
                    case 7: {
                        i.operation = "SLT";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        unsigned int rd = stoul(instruction.substr(16, 5), nullptr, 2);
                        i.r1 = rs;
                        i.r2 = rt;
                        i.r3 = rd;
                        processor->addInstruction(i);
                        break;
                    }
                    case 8: {
                        i.operation = "ADDI";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        string immediateStr = instruction.substr(16, 16) + "0000000000000000";
                        int immediate = stol(immediateStr, nullptr, 2);
                        immediate = immediate >> 16;
                        i.r1 = rs;
                        i.r2 = rt;
                        i.sValue = immediate;
                        processor->addInstruction(i);
                        break;
                    }
                    case 9: {
                        i.operation = "ANDI";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        int immediate = stol(instruction.substr(16, 16), nullptr, 2);
                        i.r1 = rs;
                        i.r2 = rt;
                        i.lValue = immediate;
                        processor->addInstruction(i);
                        break;
                    }
                    case 10: {
                        i.operation = "ORI";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        int immediate = stol(instruction.substr(16, 16), nullptr, 2);
                        i.r1 = rs;
                        i.r2 = rt;
                        i.lValue = immediate;
                        processor->addInstruction(i);
                        break;
                    }
                    case 11: {
                        i.operation = "XORI";
                        unsigned int rs = stoul(instruction.substr(6, 5), nullptr, 2);
                        unsigned int rt = stoul(instruction.substr(11, 5), nullptr, 2);
                        int immediate = stol(instruction.substr(16, 16), nullptr, 2);
                        i.r1 = rs;
                        i.r2 = rt;
                        i.lValue = immediate;
                        processor->addInstruction(i);
                        break;
                    }
                }
                continue;
            }
        }
    }
    file.close();

    processor->writeDisassemblyFile();
    processor->loadProgramData();
    processor->simulate();

    return 0;
}