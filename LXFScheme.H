
        // flux in x direction
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

        doubleField a_mix(
            0.5*(mag(U_Xpos + U_Xneg)+(a_Xpos + a_Xneg))
        );

        phix == (0.5*(
                rho_Xneg*U_Xneg + rho_Xpos*U_Xpos
                + a_mix*(rho_Xneg - rho_Xpos)
        ));

        phiUx == (0.5*(
                rho_Xneg*U_Xneg*U_Xneg+ rho_Xpos*U_Xpos*U_Xpos + p_Xpos + p_Xneg
                + a_mix*(rho_Xneg*U_Xneg - rho_Xpos*U_Xpos)
        ));

        phiEx == (0.5*(
                (rho_Xneg*(e_Xneg + 0.5*(U_Xneg*U_Xneg + V_Xneg*V_Xneg + W_Xneg*W_Xneg)) * U_Xneg + U_Xneg * p_Xneg) 
                + (rho_Xpos*(e_Xpos + 0.5*(U_Xpos*U_Xpos + V_Xpos*V_Xpos + W_Xpos*W_Xpos)) * U_Xpos + U_Xpos * p_Xpos)
                - a_mix*rho_Xpos*(e_Xpos + 0.5*(U_Xpos*U_Xpos + V_Xpos*V_Xpos + W_Xpos*W_Xpos))
                + a_mix*rho_Xneg*(e_Xneg + 0.5*(U_Xneg*U_Xneg + V_Xneg*V_Xneg + W_Xneg*W_Xneg))
        ));

        std::cout << "LXF flux calculate complete" << std::endl;