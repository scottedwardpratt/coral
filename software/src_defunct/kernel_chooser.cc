#include "kernel.h"

bool pair_match(string p1a, string p2a, string p1b, string p2b){
    string pairb(p1b+p2b);
    return (((p1a+p2a)==pairb)||((p2a+p1a)==pairb));
}

CKernel* chooseKernel(string particle1, string particle2, const parameterMap& m){
    string kdatadirname;
    bool use_cache(false);
    bool read_cache(false);
    bool hbt_only(false);

    string parsfilename = m.getS("param_filename","");
    if (m.find("use_cache")!=m.end()) {
        use_cache=m.getB("use_cache");
        read_cache = m.getB("read_cache",read_cache);
        kdatadirname = m.getS("kernel_cache","kernel_cache");
    }
    hbt_only = m.getB("hbt_only",hbt_only);
    
    CWaveFunction *wf;
    
    if (!hbt_only) {            
        if      (pair_match(particle1,particle2,"pi+","pi+")) {wf=new CWaveFunction_pipluspiplus_sqwell(parsfilename);} 
        else if (pair_match(particle1,particle2,"pi-","pi-")) {wf=new CWaveFunction_pipluspiplus_sqwell(parsfilename);} 
        else if (pair_match(particle1,particle2,"pi+","pi-")) {wf=new CWaveFunction_pipluspiminus_sqwell(parsfilename);} 
        else if (pair_match(particle1,particle2,"p","K+"))    {wf=new CWaveFunction_pkplus_sqwell(parsfilename);} 
        else if (pair_match(particle1,particle2,"p","pi+"))   {wf=new CWaveFunction_ppiplus_sqwell(parsfilename);} 
        else if (pair_match(particle1,particle2,"p","p"))     {wf=new CWaveFunction_pp_schrod(parsfilename);} 
        else if (pair_match(particle1,particle2,"K+","K+"))   {wf=new CWaveFunction_generic(parsfilename,1,KaonMass,KaonMass,1.0);} 
        else {throw MESSAGE << "Cannot compute kernel for \""<<particle1<<"\" and \""<<particle2<<"\""<<ENDM_FATAL;}

/*
        if (pair_match(particle1,particle2,"Xi","pi")) {
            kdatadirname="parameters/wfparameters_Xipi.dat";
            wf=new CWaveFunction_Xipi(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"K+","pi-")) {
            kdatadirname="parameters/wfparameters_kpluspiminus.dat";
            wf=new CWaveFunction_kpluspiminus(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"Lambda","Lambda")) {
            kdatadirname="parameters/wfparameters_lambdalambda.dat";
            wf=new CWaveFunction_lambdalambda(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"n","n")) {
            kdatadirname="parameters/wfparameters_nn.dat";
            wf=new CWaveFunction_nn(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"pi+","pi-")) {
            kdatadirname="parameters/wfparameters_pipluspiminus.dat";
            wf=new CWaveFunction_pipluspiminus(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"pi+","pi+")) {
            kdatadirname="parameters/wfparameters_pipluspiplus.dat";
            wf=new CWaveFunction_pipluspiplus(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"pi-","pi-")) {
            kdatadirname="parameters/wfparameters_piminuspiminus.dat";
            wf=new CWaveFunction_pipluspiplus(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"K+","K+")) {
            kdatadirname="parameters/wfparameters_.dat";
            wf=new CWaveFunction_pp(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"K-","K-")) {
            kdatadirname="parameters/wfparameters_.dat";
            wf=new CWaveFunction_pp(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"K+","K-")) {
            kdatadirname="parameters/wfparameters_.dat";
            wf=new CWaveFunction_pp(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"p","K+")) {
            kdatadirname="parameters/wfparameters_pkplus.dat";
            wf=new CWaveFunction_pkplus(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"p","Lambda")) {
            kdatadirname="parameters/wfparameters_plambda.dat";
            wf=new CWaveFunction_plambda(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"p","n")) {
            string kdatadirname="parameters/wfparameters_pn.dat";
            wf=new CWaveFunction_pn(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"p","p")) {
            string kdatadirname="parameters/wfparameters_pp.dat";
            wf=new CWaveFunction_pp(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"p","pi+")) {
            string kdatadirname="parameters/wfparameters_ppiplus.dat";
            wf=new CWaveFunction_ppiplus(const_cast<char*>(parsfilename.c_str()));
        } else if (pair_match(particle1,particle2,"pi0","pi0")) {
            hbt_only=true;
            kdatadirname="parameters/wfparameters_pipluspiplus.dat";
            wf=new CWaveFunction_pipluspiplus(const_cast<char*>(parsfilename.c_str()));
        } else { 
            throw MESSAGE << "Cannot compute kernel for "<<particle1<<" "<<particle2<<ENDM_FATAL;
        }
*/
    }
    
    CKernel* kernel;
    cout << endl;
    if (hbt_only)  {
        cout << "Kernel for pure HBT" <<endl;
        kernel = new CKernelExactHBT(parsfilename);
        if (!kernel) throw MESSAGE << "Kernel creation failed!"<<ENDM_FATAL;
    }
    else {
        cout << "Kernel for particle1, particle2 = "<<particle1<<", "<<particle2<<endl;
        kernel = new CKernel(parsfilename);
        if (!kernel) throw MESSAGE << "Kernel creation failed!"<<ENDM_FATAL;
        if (read_cache) {
            cout << "    Reading kernel cache... "<<endl;
            kernel->ReadData(kdatadirname);
        }
        else            kernel->Calc(wf);
        if (!read_cache && use_cache) kernel->WriteData(kdatadirname);
    }
    cout << endl;
    
    return kernel;
}
