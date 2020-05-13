% Explicit kinematic constraints of
% mg10hlDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [17x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,AE,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% jv [15x1]
%   Joint variables (rotation around z or translation in z-direction according to MDH)
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 13:01
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function jv = mg10hlDE2_kinconstr_expl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'mg10hlDE2_kinconstr_expl_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [17 1]), ...
  'mg10hlDE2_kinconstr_expl_mdh_sym_varpar: pkin has to be [17x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-11 12:59:56
% EndTime: 2020-04-11 12:59:57
% DurationCPUTime: 0.48s
% Computational Cost: add. (349->43), mult. (348->61), div. (56->7), fcn. (145->16), ass. (0->49)
t61 = -2 * pkin(8);
t41 = -pkin(8) ^ 2 + (t61 - qJ(6)) * qJ(6);
t50 = pkin(6) ^ 2 - pkin(7) ^ 2;
t60 = -t41 / 0.4e1 + t50 / 0.4e1;
t59 = -pkin(5) - pkin(4);
t58 = -pkin(5) + pkin(4);
t22 = sin(pkin(15));
t23 = cos(pkin(15));
t25 = sin(qJ(2));
t27 = cos(qJ(2));
t12 = t27 * t22 + t25 * t23;
t57 = t12 * pkin(1);
t13 = -t25 * t22 + t27 * t23;
t48 = pkin(2) * t57;
t11 = -0.2e1 * t48;
t40 = pkin(1) ^ 2;
t51 = t11 + t40;
t2 = sqrt(-((pkin(2) - t59) * (pkin(2) + t59) + t51) * ((pkin(2) - t58) * (pkin(2) + t58) + t51));
t56 = t13 * t2;
t36 = 0.1e1 / pkin(5);
t49 = pkin(2) ^ 2 + t40;
t44 = t11 + t49;
t8 = 0.1e1 / t44;
t55 = t36 * t8;
t38 = 0.1e1 / pkin(4);
t54 = t38 * t8;
t24 = qJ(6) + pkin(8);
t32 = 1 / pkin(7);
t53 = 1 / t24 * t32;
t52 = t36 * t38;
t47 = -pkin(7) - t24;
t46 = -pkin(7) + t24;
t45 = 1 / t24 ^ 2 * t32 / pkin(6);
t37 = pkin(4) ^ 2;
t43 = -t37 + t49;
t42 = (pkin(6) + t47) * (pkin(6) + t46) * (pkin(6) - t46) * (pkin(6) - t47);
t35 = pkin(5) ^ 2;
t28 = cos(pkin(16));
t26 = sin(pkin(16));
t15 = t41 + t50;
t10 = -pkin(2) * t12 + pkin(1);
t9 = -pkin(2) + t57;
t7 = sqrt(-t42);
t6 = -t35 + t37 + t44;
t5 = t11 + t35 + t43;
t4 = atan2(t7 * t53, t15 * t53);
t3 = (t15 * t60 - t42 / 0.4e1) * t45;
t1 = (-t15 / 0.4e1 + t60) * t7 * t45;
t14 = [qJ(1); qJ(2); atan2((pkin(2) * t13 * t6 + t10 * t2) * t54, (-pkin(2) * t56 + t10 * t6) * t54); -t4; atan2(t1 * t28 + t3 * t26, t1 * t26 - t3 * t28); qJ(3); qJ(4); qJ(5); t4; atan2((pkin(1) * t13 * t5 - t9 * t2) * t55, (-pkin(1) * t56 - t9 * t5) * t55); atan2(t2 * t52, (t35 - t43 + 0.2e1 * t48) * t52); atan2(t1, t3); qJ(6); 0; t61;];
jv = t14(:);
