doubleField rho_Xneg = rho.wenoRecX_neg();
doubleField rho_Xpos = rho.wenoRecX_pos();
doubleField U_Xneg = U.wenoRecX_neg();
doubleField U_Xpos = U.wenoRecX_pos();
doubleField V_Xneg = V.wenoRecX_neg();
doubleField V_Xpos = V.wenoRecX_pos();
doubleField W_Xneg = W.wenoRecX_neg();
doubleField W_Xpos = W.wenoRecX_pos();

doubleField p_Xneg = p.wenoRecX_neg();
doubleField p_Xpos = p.wenoRecX_pos();

doubleField e_Xneg = e.wenoRecX_neg();
doubleField e_Xpos = e.wenoRecX_pos();

doubleField c_Xneg = a_sonic.wenoRecX_neg();
doubleField c_Xpos = a_sonic.wenoRecX_pos();

doubleField a_Xneg(min(min(U_Xneg - c_Xneg, U_Xpos - c_Xpos), 0.*c_Xneg));
doubleField a_Xpos(max(max(U_Xneg + c_Xneg, U_Xpos + c_Xpos), 0.*c_Xneg));

phix == (
        a_Xpos*rho_Xneg*U_Xneg - a_Xneg*rho_Xpos*U_Xpos
        + a_Xpos*a_Xneg*rho_Xpos - a_Xpos*a_Xneg*rho_Xneg)
        / ( a_Xpos - a_Xneg);

phiUx == (
        a_Xpos*(rho_Xneg*U_Xneg*U_Xneg + p_Xneg)
        - a_Xneg*(rho_Xpos*U_Xpos*U_Xneg + p_Xpos)
        + a_Xpos*a_Xneg*rho_Xpos*U_Xpos 
        - a_Xpos*a_Xneg*rho_Xneg*U_Xneg)
        / (a_Xpos - a_Xneg);

phiEx == (
        a_Xpos*(rho_Xneg*(e_Xneg + 0.5*(U_Xneg*U_Xneg + V_Xneg*V_Xneg + W_Xneg*W_Xneg)) * U_Xneg + U_Xneg * p_Xneg) 
        - a_Xneg*(rho_Xpos*(e_Xpos + 0.5*(U_Xpos*U_Xpos + V_Xpos*V_Xpos + W_Xpos*W_Xpos)) * U_Xpos + U_Xpos * p_Xpos)
        + a_Xpos*a_Xneg*rho_Xpos*(e_Xpos + 0.5*(U_Xpos*U_Xpos + V_Xpos*V_Xpos + W_Xpos*W_Xpos))
        - a_Xpos*a_Xneg*rho_Xneg*(e_Xneg + 0.5*(U_Xneg*U_Xneg + V_Xneg*V_Xneg + W_Xneg*W_Xneg))
        )
        / (a_Xpos - a_Xneg);

if(rank == 0)
        std::cout << "KNP flux calculate complete" << std::endl;