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

doubleField H_Xneg(e_Xneg + p_Xneg/rho_Xneg + 0.5*(U_Xneg*U_Xneg + V_Xneg*V_Xneg + W_Xneg*W_Xneg));
doubleField H_Xpos(e_Xpos + p_Xpos/rho_Xpos + 0.5*(U_Xpos*U_Xpos + V_Xpos*V_Xpos + W_Xpos*W_Xpos));

doubleField ED_Xneg(e_Xneg + 0.5*(U_Xneg*U_Xneg + V_Xneg*V_Xneg + W_Xneg*W_Xneg));
doubleField ED_Xpos(e_Xpos + 0.5*(U_Xpos*U_Xpos + V_Xpos*V_Xpos + W_Xpos*W_Xpos));

doubleField wi_Xpos(sqrt(rho_Xpos)/(sqrt(rho_Xpos) + sqrt(rho_Xneg)));
doubleField wi_Xneg(wi_Xpos/wi_Xpos - wi_Xpos);


doubleField c_Xneg = a_sonic.wenoRecX_neg();
doubleField c_Xpos = a_sonic.wenoRecX_pos();

// doubleField aBar(0.5*(c_Xneg + c_Xpos));
// doubleField rhoBar(0.5*(rho_Xneg + rho_Xpos));

doubleField aT(wi_Xneg*c_Xneg + wi_Xpos*c_Xpos);
doubleField UT(wi_Xneg*U_Xneg + wi_Xpos*U_Xpos);

doubleField S_Xpos(min(U_Xpos - c_Xpos, UT - aT));
doubleField S_Xneg(max(U_Xneg + c_Xneg, UT + aT));

doubleField S_Star(
    (p_Xneg - p_Xpos + rho_Xpos*U_Xpos*(S_Xpos - U_Xpos) - rho_Xneg*U_Xneg*(S_Xneg - U_Xneg))
    /(rho_Xpos*(S_Xpos - U_Xpos)- rho_Xneg*(S_Xneg - U_Xneg))
);

doubleField p_Star(
    max(rho_Xpos*(U_Xpos - S_Xpos)*(U_Xpos - S_Star) + p_Xpos, 0.*p_Xpos)
);

doubleField rSDel_Xpos(rho_Xpos/rho_Xpos/(S_Xpos - S_Star));
doubleField rSDel_Xneg(rho_Xneg/rho_Xneg/(S_Xneg - S_Star));

//-- one of mass flux
doubleField phi_Xpos(rho_Xpos*U_Xpos);
doubleField phi_Xneg(rho_Xneg*U_Xneg);

doubleField phi_Xpos_Star(
    S_Star*(S_Xpos*rho_Xpos - phi_Xpos)*rSDel_Xpos
);
doubleField phi_Xneg_Star(
    S_Star*(S_Xneg*rho_Xneg - phi_Xneg)*rSDel_Xneg
);

//-- one of momentum flux
doubleField phiU_Xpos(U_Xpos*phi_Xpos + p_Xpos);
doubleField phiU_Xneg(U_Xneg*phi_Xneg + p_Xneg);

doubleField phiU_Xpos_Star(
    rSDel_Xpos*(S_Xpos - U_Xpos)*S_Star
    *(rho_Xpos*U_Xpos + (p_Star - p_Xpos)/(S_Xpos - U_Xpos))
    +p_Star
);

doubleField phiU_Xneg_Star(
    rSDel_Xneg*(S_Xneg - U_Xneg)*S_Star
    *(rho_Xneg*U_Xneg + (p_Star - p_Xneg)/(S_Xneg - U_Xneg))
    +p_Star
);

//-- one of energy flux
doubleField phiE_Xpos(H_Xpos*phi_Xpos);
doubleField phiE_Xneg(H_Xneg*phi_Xneg);

doubleField phiE_Xpos_Star(
    S_Star*(rSDel_Xpos*(S_Xpos - U_Xpos)
    *(rho_Xpos*ED_Xpos + (p_Star*S_Star - p_Xpos*U_Xpos)/(S_Xpos - U_Xpos))
    +p_Star)
);
doubleField phiE_Xneg_Star(
    S_Star*(rSDel_Xneg*(S_Xneg - U_Xneg)
    *(rho_Xneg*ED_Xneg + (p_Star*S_Star - p_Xneg*U_Xneg)/(S_Xneg - U_Xneg))
    +p_Star)
);


// doubleField a_Xneg(min(min(U_Xneg - c_Xneg, U_Xpos - c_Xpos), 0.*c_Xneg));
// doubleField a_Xpos(max(max(U_Xneg + c_Xneg, U_Xpos + c_Xpos), 0.*c_Xneg));

phix == pos0(S_Star)*(pos0(S_Xpos)*phi_Xpos
      + neg(S_Xpos)*phi_Xpos_Star)
      + neg(S_Star)*(pos0(S_Xneg)*phi_Xneg
      + neg(S_Xneg)*phi_Xneg_Star);

phiUx == pos0(S_Star)*(pos0(S_Xpos)*phiU_Xpos
      + neg(S_Xpos)*phiU_Xpos_Star)
      + neg(S_Star)*(pos0(S_Xneg)*phiU_Xneg
      + neg(S_Xneg)*phiU_Xneg_Star);

phiEx == pos0(S_Star)*(pos0(S_Xpos)*phiE_Xpos
      + neg(S_Xpos)*phiE_Xpos_Star)
      + neg(S_Star)*(pos0(S_Xneg)*phiE_Xneg
      + neg(S_Xneg)*phiE_Xneg_Star);

std::cout << "HLLC flux calculate complete" << std::endl;