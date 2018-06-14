#pragma once

#include "easylogging++.h"
// easylogging++ defines CHECK, but we want unittest version
#undef CHECK
#include "UnitTest++.h"

using namespace std;

_INITIALIZE_EASYLOGGINGPP;

/** Hard-coded paths of test data */
const string
 KPT_DUMP_DIR = "/home/hrywniak/thesis/experiments/06c_k0-with_tmat";
const string
 SI_POT_L0_PATH = "/home/hrywniak/thesis/experiments/04a_Potential_Dump/vr_jri-521_l-0.bin";
const string
 BINMAT_GPT_PATH = "/home/hrywniak/thesis/experiments/05a_K-point_Dump/gpts_001.bin";
