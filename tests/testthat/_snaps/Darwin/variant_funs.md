# run_factor

    Code
      LNR_factor[[1]]$samples$theta_lambda[, , idx]
    Output
                    F1         F2
      m      0.3821200  0.0000000
      m_lMd -0.4616134  0.5089824
      s     -0.1808061 -0.2931056
      t0     0.1766414 -0.8214164

---

    Code
      LNR_factor[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.7486287 -0.9043014
      m_lMd -0.3301907 -0.5204523
      s     -1.1308698 -0.6404446
      t0    -2.4908047 -1.7342745

---

    Code
      LNR_factor[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.3797575 -0.8663812 -0.3789995 -0.9479532 

---

    Code
      LNR_factor[[1]]$samples$theta_var[, , idx]
    Output
                      m       m_lMd           s          t0
      m      0.26335849 -0.17639169 -0.06908961  0.06749822
      m_lMd -0.17639169  0.50666538 -0.06572311 -0.49962657
      s     -0.06908961 -0.06572311  0.22048161  0.20882390
      t0     0.06749822 -0.49962657  0.20882390  0.77103439

# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9378169 -0.9489708
      m_lMd -0.3076695 -0.6012228
      s     -0.8798202 -0.4407736
      t0    -1.8683764 -1.5734354

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9534561 -0.6384648 -0.7878486 -1.7974908 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd          s        t0
      m     0.000257467 0.00000000 0.00000000 0.0000000
      m_lMd 0.000000000 0.06233054 0.00000000 0.0000000
      s     0.000000000 0.00000000 0.03554651 0.0000000
      t0    0.000000000 0.00000000 0.00000000 0.0132861

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.7791896 -0.9418038
      m_lMd -0.3359494 -0.3297219
      s     -0.9906419 -0.6841097
      t0    -2.0854534 -1.7640341

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
                m       m_lMd           s          t0 
      -1.27986980 -0.03284585 -0.59708936  0.62848861 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                   m   m_lMd          s         t0
      m     13.82729 0.00000 0.00000000 0.00000000
      m_lMd  0.00000 2.28976 0.00000000 0.00000000
      s      0.00000 0.00000 0.06131112 0.01456124
      t0     0.00000 0.00000 0.01456124 8.60142305

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8420067 -1.1170701
      m_lMd -0.2082406 -0.3745601
      s     -1.0771645 -0.3384082
      t0    -2.8310766 -1.5177944

# run_bridge

    Code
      compare(list(single = LNR_single, diag = LNR_diag, factor = LNR_factor), stage = "preburn",
      cores_for_props = 1)
    Output
               MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single -329   0  300    0  622     0        323   -23  -346 -346
      diag   -477   1 -353    1 -234     1        119  -472  -568 -592
      factor  476   0  -84    0  150     0        233  -317  -537 -550

