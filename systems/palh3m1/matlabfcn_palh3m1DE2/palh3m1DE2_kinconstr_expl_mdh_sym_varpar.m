% Explicit kinematic constraints of
% palh3m1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% jv [12x1]
%   Joint variables (rotation around z or translation in z-direction according to MDH)
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function jv = palh3m1DE2_kinconstr_expl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_kinconstr_expl_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_kinconstr_expl_mdh_sym_varpar: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 19:32:16
% EndTime: 2020-04-19 19:32:16
% DurationCPUTime: 0.38s
% Computational Cost: add. (4792->59), mult. (7192->106), div. (348->6), fcn. (4544->24), ass. (0->78)
t46 = pkin(6) ^ 2;
t51 = pkin(2) ^ 2;
t35 = sin(qJ(2));
t36 = sin(pkin(16));
t39 = cos(qJ(2));
t40 = cos(pkin(16));
t25 = t35 * t36 - t39 * t40;
t72 = pkin(5) * t25;
t59 = pkin(1) * t72;
t24 = -0.2e1 * t59;
t53 = pkin(1) ^ 2;
t62 = pkin(5) ^ 2 + t53;
t57 = t24 + t62;
t20 = -t46 + t51 + t57;
t22 = pkin(1) - t72;
t63 = t24 + t53;
t76 = -pkin(6) + pkin(2);
t77 = -pkin(6) - pkin(2);
t18 = sqrt(-((pkin(5) - t76) * (pkin(5) + t76) + t63) * ((pkin(5) - t77) * (pkin(5) + t77) + t63));
t26 = t35 * t40 + t39 * t36;
t68 = t18 * t26;
t15 = -pkin(5) * t68 + t22 * t20;
t83 = -t15 / 0.2e1;
t16 = pkin(5) * t26 * t20 + t22 * t18;
t82 = t16 / 0.2e1;
t81 = sin(pkin(17)) / 0.2e1;
t80 = sin(pkin(19)) / 0.2e1;
t79 = sin(qJ(3)) / 0.2e1;
t78 = cos(pkin(15)) / 0.2e1;
t75 = -pkin(8) - pkin(10);
t74 = -pkin(8) + pkin(10);
t38 = cos(qJ(3));
t21 = 0.1e1 / t57;
t52 = 0.1e1 / pkin(2);
t66 = t21 * t52;
t12 = (t16 * t79 + t38 * t83) * t66;
t13 = (t15 * t79 + t38 * t82) * t66;
t29 = pkin(18) + pkin(19);
t27 = sin(t29);
t28 = cos(t29);
t11 = -t28 * t12 - t27 * t13;
t73 = pkin(4) * t11;
t10 = t27 * t12 - t28 * t13;
t49 = pkin(4) ^ 2;
t60 = pkin(3) * t73;
t9 = -0.2e1 * t60;
t69 = t49 + t9;
t3 = sqrt(-((pkin(3) - t74) * (pkin(3) + t74) + t69) * ((pkin(3) - t75) * (pkin(3) + t75) + t69));
t71 = t10 * t3;
t43 = 0.1e1 / pkin(10);
t61 = pkin(3) ^ 2 + t49;
t58 = t9 + t61;
t6 = 0.1e1 / t58;
t70 = t43 * t6;
t47 = 0.1e1 / pkin(6);
t67 = t21 * t47;
t45 = 0.1e1 / pkin(8);
t65 = t43 * t45;
t64 = t47 * t52;
t56 = -t51 + t62;
t44 = pkin(8) ^ 2;
t55 = -t44 + t61;
t54 = t45 * t6 / 0.2e1;
t42 = pkin(10) ^ 2;
t37 = sin(pkin(15));
t33 = cos(pkin(19));
t31 = cos(pkin(17));
t23 = pkin(1) * t25 - pkin(5);
t19 = t24 + t46 + t56;
t17 = pkin(1) * t26 * t19 - t23 * t18;
t14 = -pkin(1) * t68 - t23 * t19;
t8 = -pkin(3) * t11 + pkin(4);
t7 = -pkin(3) + t73;
t5 = t42 + t9 + t55;
t4 = -t42 + t44 + t58;
t2 = pkin(3) * t10 * t5 + t8 * t3;
t1 = -pkin(3) * t71 + t8 * t5;
t30 = [qJ(1); qJ(2); qJ(3); atan2((t2 * t31 / 0.2e1 + t1 * t81) * t70, (-t1 * t31 / 0.2e1 + t2 * t81) * t70); qJ(4); atan2((t17 * t78 - t14 * t37 / 0.2e1) * t67, (t14 * t78 + t17 * t37 / 0.2e1) * t67); atan2((t15 * t80 + t33 * t82) * t66, (t16 * t80 + t33 * t83) * t66); atan2((pkin(4) * t10 * t4 - t7 * t3) * t54, (-pkin(4) * t71 - t7 * t4) * t54); 0; 0; atan2(t18 * t64 / 0.2e1, -(t46 - t56 + 0.2e1 * t59) * t64 / 0.2e1); atan2(t3 * t65 / 0.2e1, -(t42 - t55 + 0.2e1 * t60) * t65 / 0.2e1);];
jv = t30(:);
