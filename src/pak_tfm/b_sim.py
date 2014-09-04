'''
Created on 10/05/2012

@author: Javier Garcia Algarra

This module includes the mutualist algorithm
'''
import numpy as np
from time import time
import math
import cProfile
import sigmund_GLOBALS as sgGL
import sigmund_common as sgcom
global pathout
global invperiod
global cuentaperpl, cuentaperpol
global Logistic_abs
global model_r_alpha
global hay_bssvar
global pendiente, sd, periodo, numanyos
global count_collapse_years
global LIMIT_COLLAPSE_YEARS
global TOL_EXTINCTION

LIMIT_COLLAPSE_YEARS = 8
TOL_EXTINCTION = -0.001
invperiod = 1 / sgGL.DAYS_YEAR
count_collapse_years = 0

def signfunc(x):
    if abs(x) == x:
        return(1)
    else:
        return(-1)

def calc_coef_May(minputchar, r, K):
    return(minputchar if r != 0 else 0)

def calc_r_periodo(Ranual, invperiod):
    return(math.pow(1 + Ranual, invperiod) - 1)

def calc_r_periodo_vh(Ranual, invperiod):
    return(math.pow(1 + Ranual, invperiod) - 1)

def cuentaenlaces(mat_in):
    cuenta = 0;
    fil = len(mat_in) - 3
    col = len(mat_in[0])
    for i in range(fil):
        for j in range(col):
            if abs(mat_in[i][j]) > 0:
                cuenta += 1
    return(cuenta)

def delete_link(mat_in):
    intentos = 0
    fil = len(mat_in) - 5
    col = len(mat_in[0])
    rfil = (np.random.random_integers(0, fil - 1))
    rcol = (np.random.random_integers(0, col - 1))
    while (mat_in[rfil][rcol] == 0) and (intentos < fil * col):
        rfil = (np.random.random_integers(0, fil - 1))
        rcol = (np.random.random_integers(0, col - 1))
        intentos += 1
    mat_in[rfil][rcol] = 0
    return(mat_in, rfil, rcol)

def deletion_links_effect(k, periodoborr, minputchar_a, minputchar_b,
                          minputchar_a_mask, minputchar_b_mask, ldev_inf):
    minputchar_a, fil, col = delete_link(minputchar_a)
    minputchar_b[col][fil] = 0
    minpeq_a = minputchar_a * minputchar_a_mask
    minpeq_b = minputchar_b * minputchar_b_mask
    sgcom.inform_user(ldev_inf, \
    "Day:%d (year %d). Deleted links plant %d <-> pollinator %d"\
    % (k, k // sgGL.DAYS_YEAR , fil, col))
    return minpeq_a, minpeq_b, minputchar_a, minputchar_b,

def val_mutMay(r_species, beta, period, N1, N2, K1):
    rspneq = calc_r_periodo_vh(abs(r_species), invperiod)
    betaMay = beta * K1 / abs(r_species)
    termEq = round(betaMay * N1 * N2 / K1)
    rMay = betaMay * N2 / K1
    if (abs(termEq) > 1):
        incEq = np.random.binomial(round(abs(termEq)), -math.expm1(-rspneq)) 
    else:
        incEq = 0
    ret = [incEq, rMay]
    return(ret)

def ciclo_May(r_species, rM, inctermMay, Nindivs, K):
    rspneq = calc_r_periodo_vh(abs(r_species), invperiod)
    signosp = signfunc(r_species)
    termEq = Nindivs * (signosp - (Nindivs / K))
    rcal = r_species * ((1 - (Nindivs / K)) + rM)
    if (abs(termEq) > 1):
        incEq = np.random.binomial(round(abs(termEq)), -math.expm1(-rspneq))
    else:
        incEq = 0
    ret = [Nindivs + signfunc(termEq) * incEq + abs(inctermMay),\
           signfunc(termEq) * rcal]
    return(ret)
               
def ciclo_verhulst(rtot_species, reqsum, Nindivs, cAlpha, Alpha):
    roz = (Alpha + cAlpha * reqsum) * Nindivs
    rcal = rtot_species - roz
    rspneq = calc_r_periodo_vh(rcal, invperiod) if (rcal >= 0)\
             else calc_r_periodo_vh(-rcal, invperiod)
    incNmalth = np.random.binomial(Nindivs, -math.expm1(-rspneq))
    if (rcal > 0):
        inc_pop = Nindivs + incNmalth 
    else:
        inc_pop = Nindivs - incNmalth
    return([inc_pop, rcal])     

def ciclo_new_model(rtot_species, Nindivs, invK, L_abs):
    rcal = rtot_species
    if (rtot_species > 0):
        rcal -= (rtot_species * Nindivs * invK) 
    elif L_abs :  # |reff| correction
        rcal += (rtot_species * Nindivs * invK)      
    rspneq = calc_r_periodo_vh(rcal, invperiod) if (rcal >= 0)\
             else calc_r_periodo_vh(-rcal, invperiod)
    incNmalth = np.random.binomial(Nindivs, -math.expm1(-rspneq))
    if (rcal > 0):
        inc_pop = Nindivs + incNmalth 
    else:
        inc_pop = Nindivs - incNmalth
    return([inc_pop, rcal]) 

def ciclo_none(rtot_species, Nindivs, invK, L_abs):
    rcal = rtot_species
    rspneq = calc_r_periodo_vh(rcal, invperiod) if (rcal >= 0)\
             else calc_r_periodo_vh(-rcal, invperiod)
    incNmalth = np.random.binomial(Nindivs, -math.expm1(-rspneq))
    if (rcal > 0):
        inc_pop = Nindivs + incNmalth 
    else:
        inc_pop = Nindivs - incNmalth
    return([inc_pop, rcal])    

def calc_mutualism_params(minputchar_p, Alpha_p, Nindividuals_p, Nindividuals_q,
                          r_eqsum, term_May, rMay, k, j, n, r_p, r_muerte):
    coef_bij_matrix = float(minputchar_p[j][n] * Nindividuals_q[k][j])
    lr_eqsum = float(r_eqsum + coef_bij_matrix)
    retMay = val_mutMay(r_p[n] - r_muerte, \
            calc_coef_May(minputchar_p[j][n], r_p[n] - r_muerte, Alpha_p[n]),
                        sgGL.DAYS_YEAR, Nindividuals_p[k][n], 
                        Nindividuals_q[k][j], Alpha_p[n])
    lterm_May = term_May + retMay[0]
    lrMay = rMay + retMay[1]
    return lr_eqsum, lterm_May, lrMay

def init_lists_pop(periods, numspecies_p, minputchar_p):
    Alpha_p = []
    cAlpha_p = []
    rowNindividuals_p = []
    r_p = []
    rd_p = []
    Nindividuals_p = np.zeros([periods, numspecies_p])
    rp_eff = np.zeros([periods, numspecies_p], dtype=float)
    rp_equs = np.zeros([periods, numspecies_p], dtype=float)  
    for n in range(numspecies_p):
        rowNindividuals_p.append(int(minputchar_p[-5][n]))
        cAlpha_p.append(float(minputchar_p[-4][n]))
        Alpha_p.append(float(minputchar_p[-3][n]))
        r_p.append(minputchar_p[-2][n])
        rd_p.append(minputchar_p[-1][n])
    Nindividuals_p[0] = np.array(rowNindividuals_p)      
    return rowNindividuals_p, Alpha_p, cAlpha_p, r_p, rd_p, Nindividuals_p, \
                                rp_eff, rp_equs

def init_params_population(r_p, rd_p, n):
    r_eqsum = term_May = rMay = p_devorados = rtot_p = 0
    r_muerte = rd_p[n]
    pop_ini = 0
    return r_muerte, r_eqsum, term_May, rMay, rtot_p, p_devorados, pop_ini

def perturbation(pl_ext, n, rd_p, cuentaperp, inicioextp, periodoextp,spikep,k):
    r_m = rd_p
    cuentaper = cuentaperp
    if (k >= inicioextp) and (n in pl_ext['species']):
        posindextp = (k - inicioextp) % periodoextp
        if (posindextp >= 0) and (posindextp < spikep):
            r_m = rd_p + pl_ext['rate']
        if (posindextp == spikep - 1) and (n == pl_ext['species'][0]):
            cuentaper = cuentaper + 1
    return r_m, cuentaper

def init_forced_external_pertubations(pl_ext, pol_ext, yearperiods, 
                                      inicioextplantas, inicioextpolin, 
                                      hayextplantas, hayextpolin, 
                                      ldev_inf, lfich_inf):
    if hayextplantas:
        inicioextplantas = round(yearperiods * pl_ext['start'] * sgGL.DAYS_YEAR)
        nperpl = pl_ext['numperiod']
        periodoextpl = pl_ext['period']
        spikepl = round(periodoextpl * pl_ext['spike'])
        sgcom.inform_user(ldev_inf, "Perturbations. Plants species %s, period (years): %d, numperiods: %d, spike (fraction of period): %0.2f, rate: %.03f, start (year): %.02f" % (pl_ext['species'], periodoextpl / sgGL.DAYS_YEAR, nperpl, pl_ext['spike'], float(pl_ext['rate']), pl_ext['start']))
    else:
        inicioextplantas = nperpl = periodoextpl = spikepl = 0
    if hayextpolin:
        inicioextpolin = round(yearperiods * pol_ext['start'] * sgGL.DAYS_YEAR)
        nperpol = pol_ext['numperiod']
        periodoextpol = pol_ext['period']
        spikepol = round(periodoextpol * pol_ext['spike'])
        sgcom.inform_user(ldev_inf, "Perturbations. Pollinators species %s, period (years): %d, numperiods: %d, spike (fraction of period): %0.2f, rate: %.03f, start (year): %.02f" % (pol_ext['species'], periodoextpol / sgGL.DAYS_YEAR, nperpol, pol_ext['spike'], float(pol_ext['rate']), pol_ext['start']))
    else:
        inicioextpolin = nperpol = periodoextpol = spikepol = 0
    return nperpl, inicioextplantas, periodoextpl, spikepl, nperpol,\
           inicioextpolin, periodoextpol, spikepol

def bss_pert_sinusoidal(val):
    global periodo, sd
    modifier = 0.5 * (1.0001 + np.sin((2 * np.pi / periodo) * val))
    return np.random.normal(1, sd * modifier)

def bss_pert_lineal(val):
    modifier = 1 + pendiente * (val / numanyos)
    return np.random.normal(1, sd * modifier)

def bss_pert_no_modulation(val):
    global sd
    return(np.random.normal(1, sd))

def calcbssvarespecie(valblossomperiod, valsd, funct, nys):
    global numanyos, sd, blossomperiod
    blossomperiod = valblossomperiod
    sd = valsd
    numanyos = nys
    valdiffs = []
    for m in range (0, numanyos):
        pospl = funct(m)
        diff = max(0, 1 - ((2 / blossomperiod) * (abs(1 - pospl))))
#         if (diff < 0):
#             diff = 0
        valdiffs.append(diff)
    return valdiffs

def predators_effect(p_devorados, j, Nindividuals_p,
                     Nindividuals_c, minputchar_c, numspecies_c, n, k):
    for j in range(numspecies_c):
        if (minputchar_c[n][j] > 0):
            rceff = calc_r_periodo_vh(Nindividuals_c[k][j] * minputchar_c[n][j],
                                      invperiod)
            p_devorados = p_devorados + np.random.binomial(Nindividuals_p[k][n],
                                                        -math.expm1(-1 * rceff))
    return p_devorados, j, rceff

def populations_evolution(n, strtype, numspecies_p, algorithm, hay_foodweb, 
                          p_ext, May, haymut, rp_eff, rp_eq, 
                          minputchar_p, j, cAlpha_p, Alpha_p, r_p, rd_p, 
                          Nindividuals_p, numspecies_q, Nindividuals_q, 
                          Nindividuals_c, minputchar_c, numspecies_c, inicioext,
                          hayext, nper, periodoext, spike, k, model_r_a, 
                          ldev_inf):
    global cuentaperpl, cuentaperpol
    r_muerte, r_eqsum, term_May, rMay, rtot_p, p_devorados, pop_p = \
                                            init_params_population(r_p, rd_p, n)
    if (Nindividuals_p[k, n]):  # Much faster than (Nindividuals_p[k][n] > 0)
        # Extinciones de plantas
        if hayext:
            if (strtype == 'Plant'):
                cuentaperext = cuentaperpl  
            else: 
                cuentaperext = cuentaperpol
            if (cuentaperext < nper):
                r_muerte, cuentaperext = perturbation(p_ext, n, rd_p[n], 
                                                      cuentaperext, inicioext,
                                                      periodoext, spike, k)
                if (strtype == 'Plant'):
                    cuentaperpl = cuentaperext
                else:
                    cuentaperpol = cuentaperext
        if haymut:
            r_eqsum = np.dot(minputchar_p[0:numspecies_q, n],\
                             Nindividuals_q[k, :])
        elif (May):
            for j in range(numspecies_q):
                r_eqsum, term_May, rMay = calc_mutualism_params(minputchar_p,
                                    Alpha_p, Nindividuals_p, Nindividuals_q,
                                    r_eqsum, term_May, rMay, k, j, n, r_p,
                                    r_muerte)
        rtot_p = r_p[n] + r_eqsum - r_muerte
        
        # Efecto de los depredadores
        if hay_foodweb:
            p_devorados, j, rceff = predators_effect(p_devorados, j, 
                                                     Nindividuals_p, 
                                                     Nindividuals_c, 
                                                     minputchar_c, 
                                                     numspecies_c, n, k)
        # New algorithm
        if (model_r_a):
            retl = ciclo_verhulst(rtot_p, r_eqsum,
                                  Nindividuals_p[k, n], cAlpha_p[n], Alpha_p[n])
        elif not (May):
            retl = ciclo_new_model(rtot_p, Nindividuals_p[k, n],
                                   1 / Alpha_p[n], Logistic_abs)
        else:
            retl = ciclo_May(r_p[n] - r_muerte, rMay, term_May,
                             Nindividuals_p[k, n], Alpha_p[n])           
        pop_p = retl[0] - p_devorados
        if not(pop_p):
            sgcom.inform_user(ldev_inf,
                              "Day %d (year %d). %s species %d extincted" %\
                              (k, k // sgGL.DAYS_YEAR, strtype, n))
    Nindividuals_p[k + 1][n] = pop_p
    if (pop_p):
        rp_eff[k + 1][n], rp_eq[k + 1][n] = retl[1], rtot_p
    else:
        rp_eff[k + 1][n], rp_eq[k + 1][n] = rp_eff[k][n], rp_eq[k][n]
    rp_eff[0, ] = rp_eff[1, ]
    return 0  

def calc_compatib_plantas(numspecies, probcoinc, blossom_pert_list):
    lcomp = []
    for i in range(0, numspecies):
        if i in blossom_pert_list:
            lcomp.append(np.random.binomial(n=1, p=probcoinc))
        else:
            lcomp.append(1.0)
    return np.array(lcomp)       

def init_blossom_pertubation_params(Bssvar_period, Bssvar_sd,
                    Bssvar_species, Bssvar_modulationtype_list, numspecies_a,
                    year_periods, ldev_inf):
    # Blossom variability analysis
    global pendiente, periodo
    pBssvar_species = []
    pendiente = 0
    periodo = 0 
    hay_bssvar = (Bssvar_sd > 0.0) & (len(Bssvar_species) > 0)
    if (hay_bssvar):
        if (Bssvar_modulationtype_list[0] == 'linear'):
            strvalor = (", Slope %0.2f " % (Bssvar_modulationtype_list[1]))
        else:
            if (Bssvar_modulationtype_list[0] == 'sin'):
                strvalor = (", Period %0.2f " % (Bssvar_modulationtype_list[1]))
            else:
                strvalor = ''
        straff = (". Affected species:%s" % (Bssvar_species))
        sgcom.inform_user(ldev_inf, "Blossom variability active. Period: %0.2f (%% of year), Initial moment standard dev.: %0.04f (%% of year), Type: %s%s%s "\
                % (Bssvar_period, Bssvar_sd, Bssvar_modulationtype_list[0], 
                   strvalor, straff))
        if (str(Bssvar_species[0]).upper() == 'ALL'):
            listaspecies = list(range(numspecies_a))
        else:
            listaspecies = Bssvar_species
        bssvar_allones = np.ones((year_periods,), dtype=float)
        for i in range(numspecies_a):
            if i in listaspecies:
                if (Bssvar_modulationtype_list[0] == 'None'):
                    varspecies = calcbssvarespecie(Bssvar_period, Bssvar_sd,
                                           bss_pert_no_modulation, year_periods)
                elif (Bssvar_modulationtype_list[0] == 'linear'):
                    pendiente = float(Bssvar_modulationtype_list[1])
                    varspecies = calcbssvarespecie(Bssvar_period, Bssvar_sd,
                                           bss_pert_lineal, year_periods)
                elif (Bssvar_modulationtype_list[0] == 'sin'):
                    periodo = int(Bssvar_modulationtype_list[1])
                    varspecies = calcbssvarespecie(Bssvar_period, Bssvar_sd,
                                           bss_pert_sinusoidal, year_periods)            
                pBssvar_species.append(np.array(varspecies))
            else:
                pBssvar_species.append(bssvar_allones)
       
    return pBssvar_species, hay_bssvar, pendiente, periodo

def calc_random_blossom_effect(numspecies_a, nrows_a, ncols_a, nrows_b, ncols_b,
                               numspecies_b, plants_blossom_type, 
                               plants_blossom_prob, plants_blossom_sd, 
                               blossom_pert_list):
    if (plants_blossom_type == 'Binary'):
        lcompatibplantas = calc_compatib_plantas(numspecies_a, 
                                                 plants_blossom_prob, 
                                                 blossom_pert_list)
    else:
        lcompatibplantas = []
        for g in range(0, numspecies_a):
            if g in blossom_pert_list:
                lcompatibplantas.append(abs(np.random.normal(plants_blossom_prob, 
                                                         plants_blossom_sd, 1)))
            else:
                lcompatibplantas.append(1)
    minputchar_a_mask = np.ones([nrows_a, ncols_a])
    for i in range(0, numspecies_b):
        minputchar_a_mask[i, :] = lcompatibplantas
    minputchar_b_mask = np.ones([nrows_b, ncols_b])
    for i in range(0, numspecies_b):
        for m in range (0, numspecies_a):
            minputchar_b_mask[m, i] = lcompatibplantas[m]
    return minputchar_a_mask, minputchar_b_mask, lcompatibplantas    

def predators_population_evolution(hay_foodweb, ldev_inf, numspecies_a, 
                                   Nindividuals_a, numspecies_b, Nindividuals_b,
                                   K_c, Nindividuals_c, r_c, numspecies_c,
                                   minputchar_d, j, k):
    if hay_foodweb:
        rowNi = []
        for n in range(numspecies_c):
            coef_bij_matrix = 0
            c_devorados = 0
            r_eqsum = 0
            signo = signfunc(r_c[n])
            if (Nindividuals_c[k][n] > 0):
                for j in range(numspecies_a):
                    coef_bij_matrix = Nindividuals_a[k][j] * minputchar_d[j][n]
                    r_eqsum += coef_bij_matrix
                
                for j in range(numspecies_b):
                    coef_bij_matrix = Nindividuals_b[k][j] *\
                                      minputchar_d[j + numspecies_a][n]
                    r_eqsum += coef_bij_matrix
                
                rperiodequivalente = calc_r_periodo_vh(r_eqsum, invperiod)
                term_c = round(Nindividuals_c[k][n] *\
                               (1 - Nindividuals_c[k][n] / K_c[n]))
                if term_c > 1:
                    incPredatoria = np.random.binomial(term_c, 
                                               -math.expm1(-rperiodequivalente))
                else:
                    incPredatoria = 0
                incNmalth = np.random.binomial(Nindividuals_c[k][n], 
                        -math.expm1(-calc_r_periodo_vh(abs(r_c[n]), invperiod)))
                signo = signfunc(r_c[n])
                pop_c = Nindividuals_c[k][n] + incPredatoria + signo * incNmalth
                if (pop_c == 0):
                    sgcom.inform_user(ldev_inf, \
                            "Predator species %d extinction in day %d" % (n, k))
            else:
                pop_c = 0
            rowNi.append(pop_c)
        
        Nindividuals_c.append(rowNi)

def predators_param_init(filename, hay_foodweb, direntrada, ldev_inf, lfich_inf,
                         dt):
    K_c = []
    Nindividuals_c = []
    rowNindividuals_c = []
    r_c = []
    r_cperiod = []
    minputchar_c = []; minputchar_d = []
    numspecies_c = 0
    if hay_foodweb > 0:
        filename_c = filename + '_c.txt'
        filename_d = filename + '_d.txt'
        sgcom.inform_user(lfich_inf, "Predators matrix c:<a href='file:///" +\
                          dt + "/input/" + filename_c + "' target=_BLANK>" +\
                          filename_c + "<a><br>")
        sgcom.inform_user(lfich_inf, "Predators matrix d:<a href='file:///" +\
                          dt + "/input/" + filename_d + "' target=_BLANK>" +\
                          filename_d + "<a><br>")
        l_minputchar_c = sgcom.dlmreadlike(filename_c, direntrada)
        minputchar_c = np.array(l_minputchar_c, dtype=float)
        nrows_c = len(minputchar_c)
        ncols_c = len(minputchar_c[0])
        numspecies_c = ncols_c
        sgcom.inform_user(ldev_inf, "Predator species : %d" % numspecies_c)
        for n in range(numspecies_c):
            rowNindividuals_c.append(int(minputchar_c[nrows_c - 3][n]))
            K_c.append(int(minputchar_c[nrows_c - 2][n]))
            r_c.append(minputchar_c[nrows_c - 1][n])
            r_cperiod.append(calc_r_periodo_vh(minputchar_c[nrows_c - 1][n], 
                                               invperiod))
        Nindividuals_c.append(rowNindividuals_c)
        l_minputchar_d = sgcom.dlmreadlike(filename_d, direntrada)
        minputchar_d = np.array(l_minputchar_d, dtype=float)
        nrows_d = len(minputchar_d)
        ncols_d = len(minputchar_d[0])
        numspecies_d = ncols_d
    return Nindividuals_c, minputchar_c, numspecies_c, K_c, r_c, minputchar_d

def init_random_links_removal(eliminarenlaces, periods, ldev_inf, minputchar_a):
    periodoborr = periods
    if eliminarenlaces > 0:
        cuenta = cuentaenlaces(minputchar_a)
        hayqueborrar = math.trunc(eliminarenlaces * cuenta)
        sgcom.inform_user(ldev_inf, "Links %d. Will be deleted %d" %\
                          (cuentaenlaces(minputchar_a), hayqueborrar))
        if hayqueborrar > 0:
            periodoborr = periods // (hayqueborrar)
    return periodoborr

def init_blossom_perturbations_conditions(year_periods, blossom_pert_list, 
                                          Bssvar_period, Bssvar_sd,
                                          Bssvar_modulationtype_list,
                                          Bssvar_species, ldev_inf,
                                          numspecies_a):
    if (blossom_pert_list[0] == 'ALL'):
        bloss_species = list(range(0, numspecies_a + 1))
    else:
        bloss_species = blossom_pert_list[:]
    pBssvar_species, hay_bssvar, pendiente, periodo =\
                    init_blossom_pertubation_params(Bssvar_period, Bssvar_sd,
                                   Bssvar_species,Bssvar_modulationtype_list,
                                   numspecies_a, year_periods, ldev_inf)
    return bloss_species, hay_bssvar, pBssvar_species

def init_simulation_environment(year_periods, fichreport, algorithm, verbose):
    count_collapse_years = 0
    systemextinction = False
    periods = year_periods * sgGL.DAYS_YEAR
    May = algorithm == 'May'
    haymut = algorithm != 'NoMutualism'
    model_r_alpha = algorithm == 'Verhulst' or algorithm == 'NoMutualism'
    Logistic_abs = algorithm == 'Logistic_abs'
    ldev_inf, lfich_inf = sgcom.open_info_channels(verbose, fichreport, 'w')
    return ldev_inf, lfich_inf, periods, systemextinction, May, haymut,\
           model_r_alpha, count_collapse_years

def read_simulation_matrix(filename, dirtrabajo, direntrada, str_guild, 
                           name_guild, N0_guild, lfich_inf):
    filename_x = filename + str_guild
    dt = dirtrabajo.replace('\\', '/')
    sgcom.inform_user(lfich_inf, name_guild + " matrix: <a href='file:///" + dt +\
                      "/input/" + filename_x + "' target=_BLANK>" +\
                      filename_x + "<a>")
    l_minputchar_x = sgcom.dlmreadlike(filename_x, direntrada)
    ''' If N0_guild provided by command line'''
    if len(N0_guild) > 0:
        l_minputchar_x[-5][0] = int(N0_guild)
    minputchar_x = np.array(l_minputchar_x, dtype=float)
    try:
        nrows_x = len(minputchar_x)
    except:
        print("INPUT FILE BAD FORMAT")
    ncols_x = len(minputchar_x[0])
    numspecies_x = ncols_x
    return numspecies_x, minputchar_x, nrows_x, ncols_x


def init_blossom_perturbation_lists(plants_blossom_prob, blossom_pert_list, 
                                    numspecies_a):
    if (len(blossom_pert_list) > 0):
        lcompatibplantas = calc_compatib_plantas(numspecies_a, 
                                                 plants_blossom_prob, 
                                                 blossom_pert_list)
    else:
        lcompatibplantas = calc_compatib_plantas(numspecies_a, 
                                                 plants_blossom_prob, 
                                            [g for g in range(0, numspecies_a)])
    return lcompatibplantas

def add_report_simulation_conditions(plants_blossom_prob, plants_blossom_sd, 
                                     plants_blossom_type, blossom_pert_list,
                                     ldev_inf, numspecies_a, rowNindividuals_a,
                                     numspecies_b, rowNindividuals_b):
    sgcom.inform_user(ldev_inf, "Plant species: %d" % numspecies_a)
    sgcom.inform_user(ldev_inf, "Plant initial populations %s" %\
                      rowNindividuals_a)
    if (plants_blossom_type == 'Binary'):
        sgcom.inform_user(ldev_inf,
                "Blossom probability %s, type %s. Plant affected species:%s" %\
            (plants_blossom_prob, plants_blossom_type, str(blossom_pert_list)))
    else:
        sgcom.inform_user(ldev_inf, 
                "Blossom probability, type %s, mean %s, standard dev. %s. Plant affected species:%s" %\
                (plants_blossom_type, plants_blossom_prob, plants_blossom_sd, 
                 str(blossom_pert_list)))
    sgcom.inform_user(ldev_inf, "Pollinator species: %d" % numspecies_b)
    sgcom.inform_user(ldev_inf, 
                      "Pollinator initial populations %s" % rowNindividuals_b)

def init_external_perturbation_lists(pl_ext, pol_ext, numspecies_a, 
                                     numspecies_b):
    global cuentaperpl, cuentaperpol
    inicioextplantas = inicioextpolin = cuentaperpl = cuentaperpol = 0
    hayextplantas = len(pl_ext) > 0
    hayextpolin = len(pol_ext) > 0
    j = 0
    if hayextplantas and (pl_ext['species'][0] == 'ALL'):
        pl_ext['species'] = list(range(0, numspecies_a + 1))
    if hayextpolin and (pol_ext['species'][0] == 'ALL'):
        pol_ext['species'] = list(range(0, numspecies_b + 1))
    return inicioextplantas, inicioextpolin, hayextplantas, hayextpolin, j

def check_system_extinction(algorithm, k, lcompatibplantas, plants_blossom_prob,
                            ra_equs, rb_equs, ra_eff, rb_eff):
    global count_collapse_years
    systemextinction = False
    if (algorithm == 'Verhulst'):
        ra_xx, rb_xx = ra_eff, rb_eff
    else:
        ra_xx, rb_xx = ra_equs, rb_equs;
    # if (k>3*sgGL.DAYS_YEAR) and ((np.array(lcompatibplantas)< plants_blossom_prob).sum()==0) and \
    if (k > 3 * sgGL.DAYS_YEAR) and \
            (ra_xx[k - 1] > TOL_EXTINCTION).sum() == 0 and \
            (rb_xx[k - 1] > TOL_EXTINCTION).sum() == 0:
        count_collapse_years += 1
        if count_collapse_years >= LIMIT_COLLAPSE_YEARS:
                systemextinction = True
    else:
        count_collapse_years = max(0, count_collapse_years - 0.5)
    return systemextinction

def calc_perturbed_coeffs(k, hay_bssvar, pBssvar_species, minputchar_a, 
                          minputchar_a_mask, numspecies_a, minputchar_b, 
                          minputchar_b_mask, numspecies_b):
    minpeq_a = minputchar_a * minputchar_a_mask
    if hay_bssvar:
        for t in range(0, numspecies_b):
            for l in range (0, numspecies_a):
                minpeq_a[t, l] = minpeq_a[t, l] * pBssvar_species[l][k // sgGL.DAYS_YEAR]
    minpeq_b = minputchar_b * minputchar_b_mask
    return minpeq_a, minpeq_b

def bino_mutual(filename, year_periods, hay_foodweb, hay_superpredadores,
                data_save='', dirtrabajo='', direntrada='', dirsal='',
                eliminarenlaces=0, pl_ext=[], pol_ext=[], os='', fichreport='',
                com='', algorithm='MoMutualism', plants_blossom_prob=1.0,
                plants_blossom_sd=0.01, plants_blossom_type='Binary', 
                blossom_pert_list='', verbose=True, exit_on_extinction=False,
                N0plants='', N0pols='', release='', Bssvar_period=0.1, 
                Bssvar_sd=0.0, Bssvar_modulationtype_list=[], Bssvar_species=[]):
    global cuentaperpl, cuentaperpol
    global Logistic_abs
    global model_r_alpha
    global pendiente, blossomperiod, sd, periodo
    global count_collapse_years
    
    sgGL.ldev_inf, sgGL.lfich_inf, periods, systemextinction,\
    May, haymut, model_r_alpha, count_collapse_years = \
       init_simulation_environment(year_periods, fichreport, algorithm, verbose)    
    tinic = time()
    sgcom.start_report(sgGL.ldev_inf, filename, com, year_periods, algorithm, 
                       release)        
    numspecies_a, minputchar_a, nrows_a, ncols_a = \
                                    read_simulation_matrix(filename, dirtrabajo,
                                                           direntrada, '_a.txt',
                                                           'Plants', N0plants, 
                                                           sgGL.lfich_inf)
    rowNindividuals_a, Alpha_a, cAlpha_a, r_a, rd_a, Nindividuals_a,\
    ra_eff, ra_equs = init_lists_pop(periods, numspecies_a, minputchar_a)
    lcompatibplantas = init_blossom_perturbation_lists(plants_blossom_prob,
                                                       blossom_pert_list,
                                                       numspecies_a)
    numspecies_b, minputchar_b, nrows_b, ncols_b = \
                        read_simulation_matrix(filename, dirtrabajo, direntrada,
                                               '_b.txt', 'Pollinators ',
                                               N0pols, sgGL.lfich_inf)
    rowNindividuals_b, Alpha_b, cAlpha_b, r_b,\
    rd_b, Nindividuals_b, rb_eff, rb_equs = init_lists_pop(periods, numspecies_b,
                                                           minputchar_b)
    Nindividuals_c, minputchar_c, numspecies_c, K_c, r_c, minputchar_d = \
                                  predators_param_init(filename, hay_foodweb,
                                                      direntrada, sgGL.ldev_inf,
                                                      sgGL.lfich_inf,
                                                  dirtrabajo.replace('\\', '/'))  
    add_report_simulation_conditions(plants_blossom_prob, plants_blossom_sd,
                                     plants_blossom_type, blossom_pert_list,
                                     sgGL.ldev_inf, numspecies_a,
                                     rowNindividuals_a, numspecies_b,
                                     rowNindividuals_b)
    
    # Random links removal 
    periodoborr = init_random_links_removal(eliminarenlaces, periods, 
                                            sgGL.ldev_inf, minputchar_a)            
    # Extinction analysis. Forced death rate increases
    inicioextplantas, inicioextpolin, hayextplantas, hayextpolin, j =\
                               init_external_perturbation_lists(pl_ext, pol_ext,
                                                     numspecies_a, numspecies_b)     
    nperpl, inicioextplantas, periodoextpl, spikepl, nperpol, inicioextpolin,\
    periodoextpol, spikepol = init_forced_external_pertubations(pl_ext,
                                    pol_ext, year_periods, inicioextplantas,
                                    inicioextpolin, hayextplantas, hayextpolin,
                                    sgGL.ldev_inf, sgGL.lfich_inf)
    # Extinction analysis. Blossom pertubations   
    bloss_species, hay_bssvar, pBssvar_species = \
          init_blossom_perturbations_conditions(year_periods, blossom_pert_list,
                                                Bssvar_period, Bssvar_sd, 
                                                Bssvar_modulationtype_list,
                                                Bssvar_species, sgGL.ldev_inf,
                                                numspecies_a)
         
    for k in range (periods - 1):
        ''' The compatibilty matrixes masks are created when the year starts '''         
        if not(k % sgGL.DAYS_YEAR):  # Much faster than if ((k%sgGL.DAYS_YEAR)==0)     
            if (not(systemextinction)):
                systemextinction = check_system_extinction(algorithm, k, 
                                          lcompatibplantas, plants_blossom_prob,
                                           ra_equs, rb_equs, ra_eff, rb_eff)
                if systemextinction:
                    sgcom.inform_user(sgGL.ldev_inf,\
                         "ALARM !!!. System will collapse. Day %d (year %d)" %\
                         (k, k // sgGL.DAYS_YEAR))
                    if exit_on_extinction:
                        return(Nindividuals_a, Nindividuals_b, Nindividuals_c,
                               ra_eff, rb_eff, ra_equs, ra_equs, 
                               0, 0, 0, 0, 0, 0, systemextinction)
            minputchar_a_mask, minputchar_b_mask, lcompatibplantas =\
                      calc_random_blossom_effect(numspecies_a, nrows_a, ncols_a,
                            nrows_b, ncols_b, numspecies_b, plants_blossom_type,
                            plants_blossom_prob, plants_blossom_sd, 
                            blossom_pert_list=bloss_species[:])
            minpeq_a, minpeq_b = calc_perturbed_coeffs(k, hay_bssvar, 
                                                  pBssvar_species, minputchar_a,
                                                  minputchar_a_mask,
                                                  numspecies_a, minputchar_b,
                                                  minputchar_b_mask, 
                                                  numspecies_b)
        # Eliminacion aleatoria de enlaces
        if (eliminarenlaces > 0) & (k > 0) & (k % periodoborr == 0):
            minpeq_a, minpeq_b , minputchar_a, minputchar_b = \
                                          deletion_links_effect(k, periodoborr,
                                          minputchar_a, minputchar_b,
                                          minputchar_a_mask, minputchar_b_mask,
                                          sgGL.ldev_inf)
        [populations_evolution(n, "Plant", 
                               numspecies_a, algorithm, hay_foodweb, pl_ext, 
                               May, haymut, ra_eff, ra_equs, 
                               minpeq_a, j, cAlpha_a, Alpha_a, r_a, rd_a, 
                               Nindividuals_a, numspecies_b, Nindividuals_b, 
                               Nindividuals_c, minputchar_c, numspecies_c, 
                               inicioextplantas, hayextplantas, nperpl, 
                               periodoextpl, spikepl, k, model_r_alpha, 
                               sgGL.ldev_inf) for n in range(numspecies_a)]
        # ra_eff[0,]=ra_eff[1,]
        [populations_evolution(n, "Pollinator", 
                               numspecies_b, algorithm, hay_foodweb, pol_ext, 
                               May, haymut, rb_eff, rb_equs,
                               minpeq_b, j, cAlpha_b, Alpha_b, r_b, rd_b, 
                               Nindividuals_b, numspecies_a, Nindividuals_a,
                               Nindividuals_c, minputchar_c, numspecies_c,
                               inicioextpolin, hayextpolin, nperpol, 
                               periodoextpol, spikepol, k, model_r_alpha, 
                               sgGL.ldev_inf) for n in range(numspecies_b)]
        # rb_eff[0,]=rb_eff[1,]
        predators_population_evolution(hay_foodweb, sgGL.ldev_inf, numspecies_a,
                                       Nindividuals_a, numspecies_b, 
                                       Nindividuals_b, K_c, Nindividuals_c, r_c,
                                       numspecies_c, minputchar_d, j, k)
    maxa_individuos, maxb_individuos, max_reff, min_reff,\
    max_requs, min_requs = sgcom.find_max_values(Nindividuals_a, Nindividuals_b, 
                                           ra_eff, rb_eff, ra_equs, rb_equs)    
    tfin = time()    
    sgcom.end_report(sgGL.ldev_inf, sgGL.lfich_inf, tfin, tinic, periods, 
                     data_save, filename, algorithm, dirsal, os, Nindividuals_a,
                     ra_eff, ra_equs, Nindividuals_b, rb_eff, rb_equs, 
                     Nindividuals_c, hay_foodweb)
    return(Nindividuals_a, Nindividuals_b, Nindividuals_c, ra_eff, rb_eff,
           ra_equs, rb_equs, maxa_individuos, maxb_individuos, 
           max_reff, min_reff, max_requs, min_requs, systemextinction, 
           pBssvar_species)

            
if __name__ == '__main__':
    import doctest
    doctest.testmod()