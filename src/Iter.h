//
// Created by Anand on 23/03/15.
//

#ifndef _HYSTERESIS_ITER_H_
#define _HYSTERESIS_ITER_H_

#include <vector>
#include <tuple>
//#include "PythonCopy.h"

//struct Iterator {
//
////    std::vector<std::vector<Stack *>>::iterator list_it;
////    std::vector<Stack *>::iterator stack_it;
//
//    static std::vector<std::vector<int>> stacks;
//    static int total;
//    bool down = true;
//    int index = 0;
//
//    int * getN(int it) {
//        bool found = false;
//        int count = 0;
//
//        while (!found) {
//            if (it <= stacks[count].size() -1) {
//                found = true;
//            } else {
//                it -= stacks[count].size();
//                count++;
//            }
//        }
//
//
//        return &stacks[count][it];
//    }
//
//    int *prev(int it) {
//        if (it == 0) {
//            return &stacks[0][0];
//        }
//
//        return getN(it -1);
//    }
//
//    int *next(int it) {
//
//
//        if (it == total -1) {
//            std::vector<int> l = stacks.back();
//            int k = l.back();
//
//            return &k;
//        }
//
//        return getN(it +1);
//    }
//
//    void iter() {
//        if (index == total -1) {
//            down = false;
//        }
//
//        if (index == 0) {
//            down = true;
//        }
//
//        if (down) {
//            index++;
//        } else {
//            index--;
//        }
////        return (down) ? getN(index+1) : getN(index -1);
//    }
//
//    std::tuple<int *, int *, int *> env() {
//        return std::make_tuple(prev(index), getN(index), next(index));
//    }
//
//};

//std::vector<std::vector<int>> Iterator::stacks;


#endif //_HYSTERESIS_ITER_H_
