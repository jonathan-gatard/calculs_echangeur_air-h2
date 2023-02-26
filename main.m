%%%%%%PARAMETRES%%%%%
%Échangeur%
l_echangeur = 3.0;
d_ext_echangeur = 2.4;
d_int_echangeur = 0.8;
nombre_aller_retour = 1;
l_tube = 3.0 * nombre_aller_retour;


%Tube%
d_ext_tube = 0.010;
d_int_tube = 0.008;
section_tube = pi * ((d_int_tube / 2.0)^2.0);

%Températures%
Tinitial_air = 280.0;
Tinitial_h2 = 40.0;

%Débits%
debit_h2 = 0.51 * nombre_aller_retour;
debit_air = 8.2;


%Propriétés matériaux%
lambda_inconel = 12.0;

%Calcul%
mailles = 1000.0 * nombre_aller_retour;
deltax = l_tube / mailles;


%Tableaux%
flux_tube = zeros(1,mailles);
nb_tube_max_etage = [];
diametre_calcul = [];
Tair_tube = zeros(1,mailles);
Tair = [];
Th2 = zeros(1,mailles+1);
flux_etages = [];

%Remplissage tableaux
Tair(1)= Tinitial_air;
Th2(1)= Tinitial_h2;

%Variables boucle calcul tubes
diametre_calcul(1) = d_ext_echangeur-(d_ext_tube/2);
somme_tube=0;
nb_tube_perimetre = 0;
i = 1;

%Calcul du nombre d'étages possible
while (pi*d_int_echangeur-(d_ext_tube/2) < pi*diametre_calcul(i))
    perimetre_calcul = pi * diametre_calcul(i);
    %Calcul du nombre de tuyaux possible
    while (somme_tube <= perimetre_calcul)
        somme_tube = somme_tube + d_ext_tube + 3 * d_ext_tube;
        nb_tube_perimetre = nb_tube_perimetre + 1;
    end
    nb_tube_max_etage(i) = ceil(nb_tube_perimetre / nombre_aller_retour);
    somme_tube = 0;
    nb_tube_perimetre=0;
    diametre_calcul(i+1) = diametre_calcul(i) - 5 * d_ext_tube;
    i=i+1;
end
diametre_calcul(length(diametre_calcul)) = [];
i=1;
j=1;
k=1;
diff=1; %Précision des itérations


fprintf("Nombre d'étages : %d\n",length(nb_tube_max_etage))


while (k <= length(nb_tube_max_etage))%Sur un étage
    while (j <= mailles)%Sur un tube
        tParoiI = [];
        tParoiE = [];
        tParoiE(1) = 160.0;
        tParoiI(1) = 160.0;
        while (abs(diff) > 0.1)%Sur un deltax de tube

            %Moyenne entre la paroi extérieure et l'air
            temperatureMoyenneAir = (tParoiE(i) + Tair(k)) / 2.0;

            %Calcul des paramètre de l'air
            masseVolumiqueAir = calcul_MV_air(temperatureMoyenneAir);
            viscoDynamiqueAir = calcul_visco_air(temperatureMoyenneAir);
            surface_echangeur = (pi * diametre_calcul(k) * l_echangeur);
            vitesseAir = (debit_air) / (masseVolumiqueAir * surface_echangeur);
            reynoldsAir  =  (masseVolumiqueAir * vitesseAir * d_ext_tube) / (viscoDynamiqueAir);
            cAir  =  calcul_C_air(reynoldsAir);
            nAir = calcul_N_air(reynoldsAir);
            nusseltAir = cAir * reynoldsAir^nAir;
            lambdaAir = calcul_Lambda_air(temperatureMoyenneAir);
            cpAir = calcul_Cp_air(temperatureMoyenneAir);

            %Calcul du coefficient air
            he=((nusseltAir * lambdaAir) / (d_ext_tube));

            %Moyenne entre la paroi intérieur et l'H2
            temperatureMoyenneH2 = (tParoiI(i) + Th2(j))  /  2.0;

            %Calcul des paramètres de H2
            masseVolumiqueH2 = calcul_MV_h2(temperatureMoyenneH2);
            viscoDynamiqueH2 = calcul_visco_h2(temperatureMoyenneH2);
            vitesseH2 = (debit_h2) / (masseVolumiqueH2 * section_tube * nb_tube_max_etage(k));
            reynoldsH2 = (masseVolumiqueH2 * vitesseH2 * d_int_tube) / (viscoDynamiqueH2);
            lambdaH2 = calcul_Lambda_h2(temperatureMoyenneH2);
            cpH2 = calcul_cp_h2(temperatureMoyenneH2);
            PrandtlH2 = (viscoDynamiqueH2 * cpH2) / (lambdaH2);
            nusseltH2 = calcul_Nusselt_h2(reynoldsH2,PrandtlH2);
            
            %Calcul du coefficient H2
            hi=((nusseltH2 * lambdaH2) / (d_int_tube));


            %Surface de contact pour l'air et l'H2
            surfaceEchangeAir = (2.0 * pi * (d_ext_tube / 2.0) * deltax);
            surfaceEchangeH2 = (2.0 * pi * (d_int_tube / 2.0) * deltax);

            %Résistances équivalentes
            R1 = 1.0 / (he * surfaceEchangeAir);
            R2 = log((d_ext_tube / 2.0) / (d_int_tube / 2.0)) / (2.0 * pi * lambda_inconel * deltax);
            R3 = 1.0 / (hi * surfaceEchangeH2);
            Rth = (R1 + R2 + R3);

            %Calcul du flux réitéré (pour un deltax)
            flux  = ((Tair(k) - Th2(j)) / Rth);

            %Calcul des températures réitérées au parois (pour un deltax)
            tParoiE(i+1)=(Tair(k) - flux * R1);
            tParoiI(i+1)=(Th2(j) + flux * R3);

            %Précision du calcul
            diff=(tParoiE(i+1)-tParoiE(i));
            i=i+1;
        end
        diff=1;
        i=1;
        %Calcul sur flux sur un tube entier
        flux_tube(j) = flux;
        Th2(j+1)=((flux_tube(j)/((debit_h2/sum(nb_tube_max_etage))*cpH2)) + Th2(j));
        Tair_tube(j) = Tair(k) - (flux_tube(j)*nb_tube_max_etage(k)/(debit_air/mailles*cpAir));
        j = j + 1;          
    end
    j=1;
    if (k==1)
        plotage(Tair_tube,"Position x sur un tube","Température (K)","Variation de la témpérature de l'air le long d'un tube à d=2.4m")
    	plotage(flux_tube,"Position x sur un tube","Flux (W)","Variation du flux le long d'un tube à d=2.4m")
        fprintf('Flux de %f W sur un tube du premier étage.\n',sum(flux_tube))
    elseif (k==length(nb_tube_max_etage))
    	plotage(Tair_tube,"Position x sur un tube","Température (K)","Variation de la témpérature de l'air le long d'un tube à d=0.8m")
    	plotage(flux_tube,"Position x sur un tube","Flux (W)","Variation du flux le long d'un tube à d=0.8m")
        fprintf('Flux de %f W sur un tube du dernier étage.\n',sum(flux_tube))
    end
    
    Tair(k+1)=(Tair(k) - ((sum(flux_tube)* nb_tube_max_etage(k))/(debit_air*cpAir)));
    flux_etages(k) = sum(flux_tube) * nb_tube_max_etage(k);
    k = k + 1;
end


%Affichage températures parois

figure()
plot(tParoiI,'r')
hold on
plot(tParoiE,'g')
xlabel ('Itération')
ylabel ('Température (K)')
legend('Température paroi interieure','Température paroi extéreure')
title('Convergence de la température des parois')
hold off


plotage(Tair,"Étage","Température (K)","Variation de la température de l'air sur l'échangeur")
plotage(Th2,"Position x sur un tube","Température (K)","Variation de la température de l'H2 sur l'échangeur")
   





%Fonctions
function [] = plotage(modele,labelx,labely,titre)
	figure()
    plot(modele,'g')
    xlabel (labelx)
    ylabel (labely)
    title(titre)
    hold off
end


function [MV] = calcul_MV_air(T)
   MV=((2.7340e-5)*(T^2)) - ((1.9664e-02)*(T)) + 4.6172;
end

function [visco] = calcul_visco_air(T)
    visco =((-5.2414e-11)*(T)^2) + ((7.7422e-8)*(T)) + 6.9640e-8;
end

function [C] = calcul_C_air(R)
    if (R < 4)
        C = 0.891;
    elseif (R > 4 && R <= 40)
        C=0.821;
    elseif (R > 40 && R <= 4000)
        C=0.615;
    elseif (R > 4000 && R <= 40000)
        C=0.174;
    else
        C=0.0239;
    end
end 

function [N] = calcul_N_air(R)
    if R < 4
        N=0.330;
    elseif R > 4 && R <= 40 
        N=0.385;
    elseif R > 40 && R <= 4000
        N=0.466;
    elseif R > 4000 && R <= 40000
        N=0.618;
    else
        N=0.805;
    end
end   

function [lambda] = calcul_Lambda_air(T)
    lambda=((-2.5288e-08)*(T)^2) + ((8.9433e-05)*(T)) + 1.1043e-03; 
end     

function [Cpair] = calcul_Cp_air(T)
    if T < 233
        Cpair=0.66*T + 851.82;
    elseif T < 253
        Cpair=(-0.005*(T^2)) + 2.58 * T + 672.305;
    else
        Cpair=((-6.66667e-6)*(T^4)) + ((7.31333e-3)*(T^3)) - ((3.00779)*(T^2)) + (549.657*T) - (3.66521e4);
    end
end 

function [MV] = calcul_MV_h2(T)
   MV=2428.3*(T^(-1.124));
end

function [visco] = calcul_visco_h2(T)
    if T < 100
        visco=((1.25e-12)*(T^4)) - ((4e-10)*(T^3)) + ((4.7708e-8)*(T^2)) - (2.4957e-6)*(T) + 5.2088e-5;
    elseif T >= 100
        visco= (2.1820e-07)*(T^6.5424e-01);
    end
end
function [lambda] = calcul_Lambda_h2(T)
    if T < 40
        lambda=((-1.6417e-08)*T^4) + ((4.1830e-06)*T^3) - ((3.5711e-04)*T^2) + ((1.1094e-02)*T);
    elseif T >= 40
        lambda=2.0837e-03*T^0.79112;
    end
end
function [cpH2] = calcul_cp_h2(T)
    if (T < 50)
        cpH2=1000*(0.1785*T + 11.826);
    elseif (T >= 50)
        cpH2=1000*(((1.6125e-06)*T^4) - ((5.9195e-04)*T^3) + ((8.2675E-02)*T^2) - (5.2271*T) + 139.34);
    elseif (T >= 100)
        cpH2=1000*(((-3.2125e-11)*T^4) + ((8.5442e-08)*T^3) - ((7.8702e-05)*T^2) + ((3.0481e-02)*T) + 10.312);
    end
end
function [Nu] = calcul_Nusselt_h2(R,P)
    if(R<2500)
        Nu=4.363;
    elseif(R>=2500)
        Nu=0.023*(R^0.8)*(P^(1/3));
    end
end

