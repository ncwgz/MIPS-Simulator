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
#include <queue>
#include <algorithm>

using namespace std;

struct Instruction {
    string      assembly;
    int         category;
    int         opcode;
    string      operation;
    unsigned    int r1 = -1;
    unsigned    int r2 = -1;
    unsigned    int r3 = -1;
    long        lValue;
    int         sValue;
    int         index;
};

bool operator < (const Instruction &i1, const Instruction &i2) {
    return (i1.index < i2.index);
}

bool operator == (const Instruction &i1, const Instruction &i2) {
    return (i1.index == i2.index);
}

struct WriteData {
    Instruction i;
    int p_order;
    int r;
    long value;
};

class Processor {
private:
    unsigned int            PC;
    long*                   registers;
    vector<int>             r_ready[32];
    vector<int>             w_ready[32];
    map<int, long>*         memory;
    map<int, Instruction>*  instructions;
    vector<int>*            altered_memory;
    int                     cycle_counter;

    vector<pair<int, Instruction>>*    pre_issue;  // 4 entries
    vector<pair<int, Instruction>>*    pre_alu_1;  // 2 entries
    vector<pair<int, Instruction>>*    pre_alu_2;  // 2 entries
    vector<WriteData>*                 post_alu_2; // 1 entry
    vector<pair<int, Instruction>>*    pre_mem;    // 1 entry
    vector<WriteData>*                 post_mem;   // 1 entry

    bool isToBreak = false;
    Instruction* waiting_instruction = nullptr;
    Instruction* executed_instruction = nullptr;
    int fetch_counter = 100;

    void lock (Instruction i, int p_order) {
        if (i.category == 1) {
            // mem[r1 + value] <- r2
            if (i.opcode == 6) {
                w_ready[i.r1].push_back(p_order);
                w_ready[i.r2].push_back(p_order);
            }
            // r2 = r1 op value
            else {
                r_ready[i.r2].push_back(p_order);
                w_ready[i.r2].push_back(p_order);

                w_ready[i.r1].push_back(p_order);
            }
        }
        if (i.category == 3) {
            // r3 <- r1 op r2
            if (i.opcode < 8) {
                r_ready[i.r3].push_back(p_order);
                w_ready[i.r3].push_back(p_order);

                w_ready[i.r1].push_back(p_order);
                if (i.r1 != i.r2)
                    w_ready[i.r2].push_back(p_order);
            }
            // r2 <- r1 op value
            else {
                r_ready[i.r2].push_back(p_order);
                w_ready[i.r2].push_back(p_order);

                w_ready[i.r1].push_back(p_order);
            }
        }
    }

    // used in unlock
    void erase_lock (bool is_read, int reg, int p_order) {
        if (is_read) {
            r_ready[reg].erase(remove(r_ready[reg].begin(), r_ready[reg].end(),
                            p_order), r_ready[reg].end());
        } else {
            w_ready[reg].erase(remove(w_ready[reg].begin(), w_ready[reg].end(),
                            p_order), w_ready[reg].end());
        }
    }

    void unlock (Instruction i, int p_order) {
        if (i.category == 1) {
            // mem[r1 + value] <- r2
            if (i.opcode == 6) {
                erase_lock(false, i.r1, p_order);
                erase_lock(false, i.r2, p_order);
            }
            // r2 = r1 op value
            else {
                erase_lock(true, i.r2, p_order);
                erase_lock(false, i.r2, p_order);

                erase_lock(false, i.r1, p_order);
            }
        }
        if (i.category == 3) {
            // r3 <- r1 op r2
            if (i.opcode < 8) {
                erase_lock(true, i.r3, p_order);
                erase_lock(false, i.r3, p_order);

                erase_lock(false, i.r1, p_order);
                erase_lock(false, i.r2, p_order);
            }
            // r2 <- r1 op value
            else {
                erase_lock(true, i.r2, p_order);
                erase_lock(false, i.r2, p_order);

                erase_lock(false, i.r1, p_order);
            }
        }
    }

    // return true if locked
    bool check_lock (bool is_read, int reg, int p_order) {
        if (is_read) {
            if (!r_ready[reg].empty()) {
                if (r_ready[reg].front() < p_order) {
                    return true;
                }
            }
        } else {
            if (!w_ready[reg].empty()) {
                if (w_ready[reg].front() < p_order) {
                    return true;
                }
            }
        }
        return false;
    }

    // FETCH
    vector<pair<int, Instruction>>* fetch() {
        auto to_fetch = new vector<pair<int, Instruction>>();

        this->executed_instruction = nullptr;
        for (int i = 0; i < 2; ++i) {
            auto inst = instructions->find(PC)->second;
            // pre_issue_queue is full
            if (pre_issue->size() + to_fetch->size() == 4) break;
            // Non-branch operations
            if (inst.category == 3 || (inst.category == 1 && inst.opcode > 5 && inst.opcode != 11)) {
                to_fetch->push_back(make_pair(fetch_counter, inst));
                lock (inst, fetch_counter);
                fetch_counter ++;
                PC += 4;
                continue;
            }
            // NOP
            if (inst.opcode == 11) {
                this->executed_instruction = &instructions->find(PC)->second;
                PC += 4;
                continue;
            }
            // BREAK
            if (inst.opcode == 5) {
                this->isToBreak = true;
                Instruction ii = inst;
                this->executed_instruction = &instructions->find(PC)->second;
                break;
            }
            // J
            if (inst.opcode == 0) {
                this->executed_instruction = &instructions->find(PC)->second;
                PC = inst.lValue;
                break;
            }
            // JR, BLTZ, BLGZ
            if (inst.opcode == 1 || inst.opcode == 3 || inst.opcode == 4) {
                if (!check_lock(true, inst.r1, fetch_counter)) {
                    switch (inst.opcode) {
                        // JR
                        case 1: {
                            this->executed_instruction = &instructions->find(PC)->second;
                            PC = registers[inst.r1];
                            this->waiting_instruction = nullptr;
                            break;
                        }
                            // BLTZ
                        case 3: {
                            if (registers[inst.r1] < 0) {
                                this->executed_instruction = &instructions->find(PC)->second;
                                PC += inst.lValue;
                                this->waiting_instruction = nullptr;
                            }
                            break;
                        }
                            // BGTZ
                        case 4: {
                            if (registers[inst.r2] > 0) {
                                this->executed_instruction = &instructions->find(PC)->second;
                                PC += inst.lValue;
                                this->waiting_instruction = nullptr;
                            }
                            break;
                        }
                    }
                    break;
                } else {
                    this->waiting_instruction = &instructions->find(PC)->second;
                    break;
                }
            }
            // BEQ
            if (inst.opcode == 2) {
                if ((!check_lock(true, inst.r1, fetch_counter)) && (!check_lock(true,inst.r2, fetch_counter))) {
                    this->executed_instruction = &instructions->find(PC)->second;
                    this->waiting_instruction = nullptr;
                    if (registers[inst.r1] == registers[inst.r2]) {
                        PC += inst.lValue;
                    } else {
                        PC += 4;
                    }
                    break;
                } else {
                    this->waiting_instruction = &instructions->find(PC)->second;
                    break;
                }
            }
        }
        return to_fetch;
    }
    void ffetch(vector<pair<int, Instruction>> to_fetch) {
        for (auto it : to_fetch) {
            pre_issue->push_back(it);
        }
    }

    // ISSUE
    vector<pair<int, Instruction>>* issue() {
        bool is_alu1_allocated = false;
        bool is_alu2_allocated = false;
        auto *to_issue = new vector<pair<int, Instruction>>();

        auto it = pre_issue->begin();
        while (it != pre_issue->end()) {

            // ALU1 & ALU2 both allocated
            if (is_alu1_allocated && is_alu2_allocated) {
                break;
            }

            Instruction i = it->second;
            int o = it->first;

            // ======= check registers =======
            // SW: read r1, r2
            if (i.category == 1 && i.opcode == 6) {
                if (check_lock(true, i.r1, o) || check_lock(true, i.r2, o)) {
                    it++;
                    continue;
                }
            }
            // write r3, read r1, r2
            else if (i.category == 3 && i.opcode < 8) {
                if (check_lock(false, i.r3, o) || check_lock(true, i.r1, o) || check_lock(true, i.r2, o)) {
                    it++;
                    continue;
                }
            }
            // write r2, read r1
            else {
                if (check_lock(false, i.r2, o) || check_lock(true, i.r1, o)) {
                    it++;
                    continue;
                }
            }
            // ============ done =============

            // ========= check stores ========
            // LW/SW: No SW before it
            if (i.category == 1 && i.opcode < 8) {
                for (auto itt : *pre_issue) {
                    if (itt.second.category == 1 && itt.second.opcode == 6) {
                        if (itt.first < it->first) {
                            it ++;
                            continue;
                        }
                    }
                }
            }
            // ============ done =============

            // ALU2
            if (i.category == 3 || (i.category == 1 && i.opcode > 7)) {
                if (!is_alu2_allocated) {
                    to_issue->push_back(*it);
                    is_alu2_allocated = true;
                }
                it ++;
                continue;
            }
            // ALU1
            else {
                if (!is_alu1_allocated) {
                    to_issue->push_back(*it);
                    is_alu1_allocated = true;
                }
                it ++;
                continue;
            }
        }

        for (auto it : *to_issue) {
            pre_issue->erase(
                    remove(pre_issue->begin(), pre_issue->end(), it),
                    pre_issue->end());
        }
        return to_issue;
    }
    void iissue(vector<pair<int, Instruction>> to_issue) {
        for (auto it : to_issue) {
            // ALU2
            if (it.second.category == 3 || (it.second.category == 1 && it.second.opcode > 7)) {
                pre_alu_2->push_back(it);
            }
            // ALU1
            else {
                pre_alu_1->push_back(it);
            }
        }
    }

    // ALU1
    pair<int, Instruction>* alu1() {
        pair<int, Instruction>* to_mem = nullptr;
        if (!pre_alu_1->empty()) {
            auto *p = new pair<int, Instruction>();
            *p = pre_alu_1->front();
            to_mem = p;
            pre_alu_1->erase(pre_alu_1->begin());
        }
        return to_mem;
    }
    void aalu1(pair<int, Instruction>* to_mem) {
        if (to_mem) {
            pre_mem->push_back(*to_mem);
        }
    }

    // ALU2
    WriteData* alu2() {
        WriteData* to_write = nullptr;
        if (!pre_alu_2->empty()) {
            auto it = pre_alu_2->front();
            Instruction i = it.second;
            to_write = new WriteData;
            to_write->i = i;
            to_write->p_order = it.first;

            if (i.category == 1) {
                to_write->r = i.r2;
                switch (i.opcode) {
                    case 8: { // SLL
                        to_write->value = registers[i.r1] << i.lValue; break;
                    }
                    case 9: { // SRL
                        to_write->value = ((unsigned int)registers[i.r1]) >> i.lValue; break;
                    }
                    case 10: { // SRA
                        to_write->value = registers[i.r1] >> i.lValue; break;
                    }
                }
            }
            else if (i.category == 3 && i.opcode < 8) {
                to_write->r = i.r3;
                switch (i.opcode) {
                    case 0: { // ADD
                        to_write->value = registers[i.r1] + registers[i.r2]; break;
                    }
                    case 1: { // SUB
                        to_write->value = registers[i.r1] - registers[i.r2]; break;
                    }
                    case 2: { // MUL
                        to_write->value = registers[i.r1] * registers[i.r2]; break;
                    }
                    case 3: { // AND
                        to_write->value = registers[i.r1] & registers[i.r2]; break;
                    }
                    case 4: { // OR
                        to_write->value = registers[i.r1] | registers[i.r2]; break;
                    }
                    case 5: { // XOR
                        to_write->value = registers[i.r1] ^ registers[i.r2]; break;
                    }
                    case 6: { // NOR
                        to_write->value = ~(registers[i.r1] | registers[i.r2]); break;
                    }
                    case 7: { // SLT
                        if (registers[i.r1] < registers[i.r2]) {
                            to_write->value = 1;
                        } else {
                            to_write->value = 0;
                        }
                        break;
                    }
                }
            }
            else if (i.category == 3 && i.opcode > 7) {
                to_write->r = i.r2;
                switch (i.opcode) {
                    case 8: { // ADDI
                        to_write->value = registers[i.r1] + i.sValue; break;
                    }
                    case 9: { // ANDI
                        to_write->value = registers[i.r1] & i.lValue; break;
                    }
                    case 10: { // ORI
                        to_write->value = registers[i.r1] | i.lValue; break;
                    }
                    case 11: { // XORI
                        to_write->value = registers[i.r1] ^ i.lValue; break;
                    }
                }
            }

            pre_alu_2->erase(pre_alu_2->begin());
        }
        return to_write;
    }
    void aalu2(WriteData *to_write) {
        if (to_write) {
            post_alu_2->push_back(*to_write);
        }
    }

    // MEM
    WriteData* mem() {
        WriteData* to_write = nullptr;
        if (!pre_mem->empty()) {
            Instruction i = (*pre_mem)[0].second;
            // SW: Do now
            if (i.opcode == 6) {
                memory->find(registers[i.r1] + i.lValue)->second = registers[i.r2];
            }
            // LW: Generate WriteData
            else {
                to_write = new WriteData;
                to_write->i = i;
                to_write->p_order = (*pre_mem)[0].first;
                to_write->r = i.r2;
                to_write->value = memory->find(registers[i.r1] + i.lValue)->second;
            }
            pre_mem->erase(pre_mem->begin());
        }
        return to_write;
    }
    void mmem(WriteData* to_write) {
        if (to_write) {
            post_mem->push_back(*to_write);
        }
    }

    // WB
    void writeback() {
        if (!post_mem->empty()) {
            WriteData data = post_mem->front();
            registers[data.r] = data.value;
            post_mem->erase(post_mem->begin());
            unlock(data.i, data.p_order);
        }
        if (!post_alu_2->empty()) {
            WriteData data = post_alu_2->front();
            registers[data.r] = data.value;
            post_alu_2->erase(post_alu_2->begin());
            unlock(data.i, data.p_order);
        }
    }

public:
    Processor () {
        this->PC = 256;
        this->registers = new long[32];
        this->memory = new map<int, long>();
        this->instructions = new map<int, Instruction>();
        this->altered_memory = new vector<int>();

        for (int i = 0; i < 32; ++i) {
            registers[i] = 0;
        }

        this->cycle_counter = 1;

        this->pre_issue = new vector<pair<int, Instruction>>();
        this->pre_alu_1 = new vector<pair<int, Instruction>>();
        this->pre_alu_2 = new vector<pair<int, Instruction>>();
        this->post_alu_2 = new vector<WriteData>();
        this->pre_mem = new vector<pair<int, Instruction>>();
        this->post_mem = new vector<WriteData>();
    }

    void addInstruction(Instruction i) {
        this->instructions->insert(pair<int, Instruction>(this->PC, i));
        PC += 4;
    }

    string getAssembly(Instruction i, bool isIndexNeeded = true) {
        stringstream ss;
        if (isIndexNeeded) {
            ss << i.index << "\t";
        }
        if (i.category == 0) {
            ss << i.sValue;
        }
        else if (i.category == 1) {
            ss << i.operation << " ";
            switch (i.opcode) {
                case 0: { // J
                    ss << "#" << i.lValue;
                    break;
                }
                case 1: { // JR
                    ss << "R" << i.r1;
                    break;
                }
                case 2: { // BEQ
                    ss << "R" << i.r1 << ", ";
                    ss << "R" << i.r2 << ", ";
                    ss << "#" << i.lValue;
                    break;
                }
                case 3: { // BLGZ
                    // the same to case 4(BGTZ)
                }
                case 4: { // BGTZ
                    ss << "R" << i.r1 << ", ";
                    ss << "#" << i.lValue;
                    break;
                }
                case 5: { // BREAK
                    break;
                }
                case 6: { // SW
                    // the same to case 7(LW)
                }
                case 7: { // LW
                    ss << "R" << i.r2 << ", ";
                    ss << "" << i.lValue;
                    ss << "(R" << i.r1 << ")";
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
                    ss << "#" << i.lValue;
                    break;
                }
                case 11: { // NOP
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
                    ss << "R" << i.r2;
                    break;
                }
                case 8: { // ADDI
                    ss << "R" << i.r2 << ", ";
                    ss << "R" << i.r1 << ", ";
                    ss << "#" << i.sValue;
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
                    ss << "#" << i.lValue;
                    break;
                }
            }
        }
        if (isIndexNeeded) {
            ss << endl;
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

    void writeSimulationFile(ofstream& file) {
        file << "--------------------" << endl;
        file << "Cycle " << cycle_counter << ":" << endl;
        file << endl;

        file << "IF Unit:" << endl;
        file << "\tWaiting Instruction:";
        if (waiting_instruction)
            file << " [" << getAssembly(*waiting_instruction, false) << "]";
        file << endl;
        file << "\tExecuted Instruction:";
        if (executed_instruction)
            file << " [" << getAssembly(*executed_instruction, false) << "]";
        file << endl;

        file << "Pre-Issue Queue:" << endl;
        for (int i = 0; i < 4; ++i) {
            file << "\tEntry " << i << ":";
            if (pre_issue->size() > i) {
                file << " [" << getAssembly((*pre_issue)[i].second, false) << "]";
            }
            file << endl;
        }

        file << "Pre-ALU1 Queue:" << endl;
        for (int i = 0; i < 2; ++i) {
            file << "\tEntry " << i << ":";
            if (pre_alu_1->size() > i) {
                    file << " [" << getAssembly((*pre_alu_1)[i].second, false) << "]";
            }
            file << endl;
        }

        file << "Pre-MEM Queue:";
        if (!pre_mem->empty()) {
            file << " [" << getAssembly(pre_mem->front().second, false) << "]";
        }
        file << endl;

        file << "Post-MEM Queue:";
        if (!post_mem->empty()) {
            file << " [" << getAssembly(post_mem->front().i, false) << "]";
        }
        file << endl;

        file << "Pre-ALU2 Queue:" << endl;
        for (int i = 0; i < 2; ++i) {
            file << "\tEntry " << i << ":";
            if (pre_alu_2->size() > i) {
                file << " [" << getAssembly((*pre_alu_2)[i].second, false) << "]";
            }
            file << endl;
        }

        file << "Post-ALU2 Queue:";
        if (!post_alu_2->empty()) {
            file << " [" << getAssembly(post_alu_2->front().i, false) << "]";
        }
        file << endl;

        file << endl;
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
        this->cycle_counter = 0;
        this->altered_memory->clear();
        bool isToBreak = false;
        this->PC = 256;
        ofstream file("simulation.txt", ios::out);
        while (true) {
            cycle_counter ++;
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

            writeSimulationFile(file);
            if (isToBreak) break;
        }
        file << endl;
        file.close();
    }

    void simulate_with_pipeline() {
        this->cycle_counter = 0;
        this->altered_memory->clear();
        this->PC = 256;
        ofstream file("simulation.txt", ios::out);
        while (!isToBreak) {
            cycle_counter ++;
//            cout << cycle_counter << ", " << PC << endl;

            // At the beginning of the cycle:
            auto to_fetch = fetch();
            auto to_issue = issue();
            auto to_mem = alu1();
            auto to_write1 = mem();
            auto to_write2 = alu2();
            writeback();

            // At the end of the cycle:
            ffetch(*to_fetch);
            iissue(*to_issue);
            aalu1(to_mem);
            mmem(to_write1);
            aalu2(to_write2);

            writeSimulationFile(file);

            cout << "Cycle " << cycle_counter << ": ============" << endl;
            cout << "R" << endl;
            for (int i = 0; i < 32; ++i) {
                auto it = r_ready[i];
                if (!it.empty()) {
                    cout << "\t#";
                    cout << i << ":";
                    for (auto itt : it) {
                        cout << itt << "\t";
                    }
                    cout << endl;
                }
            }
            cout << "W" << endl;
            for (int i = 0; i < 32; ++i) {
                auto it = w_ready[i];
                if (!it.empty()) {
                    cout << "\t#";
                    cout << i << ":";
                    for (auto itt : it) {
                        cout << itt << "\t";
                    }
                    cout << endl;
                }
            }
            cout << endl;
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
    processor->simulate_with_pipeline();

    return 0;
}