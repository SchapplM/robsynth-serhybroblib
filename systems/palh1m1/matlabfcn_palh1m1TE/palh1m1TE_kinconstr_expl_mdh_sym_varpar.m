% Explicit kinematic constraints of
% palh1m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% jv [16x1]
%   Joint variables (rotation around z or translation in z-direction according to MDH)
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function jv = palh1m1TE_kinconstr_expl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_kinconstr_expl_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_kinconstr_expl_mdh_sym_varpar: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:08:06
% EndTime: 2020-04-12 20:08:07
% DurationCPUTime: 0.48s
% Computational Cost: add. (4920->84), mult. (7383->144), div. (360->9), fcn. (4665->31), ass. (0->107)
t59 = pkin(8) ^ 2;
t65 = pkin(3) ^ 2;
t46 = sin(qJ(2));
t47 = sin(pkin(19));
t50 = cos(qJ(2));
t51 = cos(pkin(19));
t34 = t46 * t51 - t50 * t47;
t97 = pkin(7) * t34;
t79 = pkin(1) * t97;
t31 = -0.2e1 * t79;
t69 = pkin(1) ^ 2;
t84 = pkin(7) ^ 2 + t69;
t76 = t31 + t84;
t23 = -t59 + t65 + t76;
t29 = pkin(1) - t97;
t102 = -pkin(8) + pkin(3);
t103 = -pkin(8) - pkin(3);
t85 = t31 + t69;
t19 = sqrt(-((pkin(7) - t102) * (pkin(7) + t102) + t85) * ((pkin(7) - t103) * (pkin(7) + t103) + t85));
t35 = t46 * t47 + t50 * t51;
t92 = t19 * t35;
t15 = -pkin(7) * t92 + t29 * t23;
t112 = -t15 / 0.2e1;
t16 = pkin(7) * t35 * t23 + t29 * t19;
t111 = t16 / 0.2e1;
t110 = sin(pkin(23)) / 0.2e1;
t109 = sin(pkin(21)) / 0.2e1;
t49 = cos(qJ(3));
t108 = -t49 / 0.2e1;
t107 = cos(pkin(18)) / 0.2e1;
t106 = 0.1e1 / pkin(2) / 0.2e1;
t105 = -pkin(2) - pkin(13);
t104 = -pkin(2) + pkin(13);
t101 = -pkin(9) - pkin(11);
t100 = -pkin(9) + pkin(11);
t45 = sin(qJ(3));
t25 = 0.1e1 / t76;
t66 = 0.1e1 / pkin(3);
t89 = t25 * t66;
t12 = (t16 * t108 + t45 * t112) * t89;
t13 = (t15 * t108 + t45 * t111) * t89;
t38 = pkin(23) + pkin(22);
t36 = sin(t38);
t37 = cos(t38);
t10 = t37 * t12 + t36 * t13;
t99 = pkin(5) * t10;
t41 = sin(pkin(20));
t44 = cos(pkin(20));
t32 = -t49 * t41 - t45 * t44;
t98 = pkin(6) * t32;
t11 = -t36 * t12 + t37 * t13;
t63 = pkin(5) ^ 2;
t81 = pkin(4) * t99;
t9 = -0.2e1 * t81;
t94 = t63 + t9;
t3 = sqrt(-((pkin(4) - t100) * (pkin(4) + t100) + t94) * ((pkin(4) - t101) * (pkin(4) + t101) + t94));
t96 = t11 * t3;
t56 = 0.1e1 / pkin(11);
t82 = pkin(4) ^ 2 + t63;
t78 = t9 + t82;
t6 = 0.1e1 / t78;
t95 = t56 * t6;
t80 = pkin(1) * t98;
t28 = -0.2e1 * t80;
t62 = pkin(6) ^ 2;
t86 = t28 + t62;
t18 = sqrt(-((pkin(1) - t104) * (pkin(1) + t104) + t86) * ((pkin(1) - t105) * (pkin(1) + t105) + t86));
t33 = t45 * t41 - t49 * t44;
t93 = t18 * t33;
t83 = t62 + t69;
t77 = t28 + t83;
t24 = 0.1e1 / t77;
t54 = 0.1e1 / pkin(13);
t91 = t24 * t54;
t60 = 0.1e1 / pkin(8);
t90 = t25 * t60;
t58 = 0.1e1 / pkin(9);
t88 = t56 * t58;
t87 = t60 * t66;
t75 = -t65 + t84;
t67 = pkin(2) ^ 2;
t74 = -t67 + t83;
t57 = pkin(9) ^ 2;
t73 = -t57 + t82;
t72 = t58 * t6 / 0.2e1;
t71 = t24 * t106;
t70 = t54 * t106;
t55 = pkin(11) ^ 2;
t53 = pkin(13) ^ 2;
t48 = sin(pkin(18));
t43 = cos(pkin(21));
t42 = cos(pkin(23));
t30 = pkin(1) * t34 - pkin(7);
t27 = -pkin(1) * t32 + pkin(6);
t26 = -pkin(1) + t98;
t22 = t31 + t59 + t75;
t21 = t28 + t53 + t74;
t20 = -t53 + t67 + t77;
t17 = pkin(1) * t35 * t22 - t30 * t19;
t14 = -pkin(1) * t92 - t30 * t22;
t8 = -pkin(4) * t10 + pkin(5);
t7 = -pkin(4) + t99;
t5 = t55 + t9 + t73;
t4 = -t55 + t57 + t78;
t2 = pkin(4) * t11 * t5 + t8 * t3;
t1 = -pkin(4) * t96 + t8 * t5;
t39 = [qJ(1); qJ(2); qJ(3); atan2((t1 * t109 + t2 * t43 / 0.2e1) * t95, (-t1 * t43 / 0.2e1 + t2 * t109) * t95); qJ(4); atan2((t17 * t107 + t14 * t48 / 0.2e1) * t90, (t14 * t107 - t48 * t17 / 0.2e1) * t90); atan2((t15 * t110 + t42 * t111) * t89, (t16 * t110 + t42 * t112) * t89); atan2((pkin(6) * t33 * t20 - t26 * t18) * t71, (-pkin(6) * t93 - t26 * t20) * t71); atan2(t18 * t70, (t53 - t74 + 0.2e1 * t80) * t70); atan2((pkin(5) * t11 * t4 - t7 * t3) * t72, (-pkin(5) * t96 - t7 * t4) * t72); atan2(t19 * t87 / 0.2e1, -(t59 - t75 + 0.2e1 * t79) * t87 / 0.2e1); atan2((pkin(1) * t33 * t21 + t27 * t18) * t91 / 0.2e1, -(-pkin(1) * t93 + t27 * t21) * t91 / 0.2e1); atan2(t3 * t88 / 0.2e1, -(t55 - t73 + 0.2e1 * t81) * t88 / 0.2e1); 0; 0; 0;];
jv = t39(:);
