% Calculate joint inertia matrix for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2DE2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_inertiaJ_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_inertiaJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE2_inertiaJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE2_inertiaJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:58:21
% EndTime: 2020-05-02 20:58:22
% DurationCPUTime: 1.08s
% Computational Cost: add. (1531->163), mult. (2651->252), div. (0->0), fcn. (3524->22), ass. (0->81)
t72 = cos(qJ(4));
t59 = t72 ^ 2;
t67 = sin(qJ(4));
t95 = t67 ^ 2 + t59;
t68 = sin(qJ(3));
t73 = cos(qJ(3));
t112 = ((mrSges(11,2) + mrSges(4,2)) * t73 + (mrSges(11,1) + mrSges(4,1)) * t68) * pkin(1);
t65 = cos(pkin(20));
t70 = sin(pkin(18));
t75 = cos(pkin(18));
t94 = sin(pkin(20));
t42 = t75 * t65 + t70 * t94;
t90 = t65 * t70 - t75 * t94;
t27 = -t68 * t42 + t90 * t73;
t31 = t42 * t73 + t68 * t90;
t69 = sin(qJ(2));
t74 = cos(qJ(2));
t15 = t27 * t74 - t69 * t31;
t16 = t69 * t27 + t31 * t74;
t93 = pkin(22) + pkin(21);
t54 = cos(t93);
t86 = sin(t93);
t111 = t15 * t54 - t86 * t16;
t104 = pkin(1) * t73;
t9 = -t15 * t86 - t16 * t54;
t92 = pkin(1) * t68 + pkin(5);
t2 = t111 * t104 - t9 * t92;
t110 = t2 ^ 2;
t109 = -2 * pkin(15);
t63 = sin(pkin(19));
t66 = cos(pkin(19));
t38 = t63 * t73 + t66 * t68;
t39 = -t63 * t68 + t66 * t73;
t28 = -t38 * t69 + t39 * t74;
t22 = -pkin(2) * t28 - pkin(15);
t108 = 0.2e1 * t22;
t106 = t2 * t9;
t4 = -t9 * t104 - t111 * t92;
t105 = t4 * t111;
t103 = pkin(2) * mrSges(10,3);
t78 = pkin(5) ^ 2;
t102 = t111 ^ 2 * t78;
t20 = t42 * t86 - t90 * t54;
t101 = mrSges(6,3) * t20;
t100 = t67 * Ifges(6,6);
t97 = t72 * Ifges(6,5);
t21 = -t42 * t54 - t90 * t86;
t96 = -Ifges(6,3) * t21 + t20 * t97;
t53 = t69 * pkin(1) - pkin(15);
t85 = mrSges(6,1) * t67 + mrSges(6,2) * t72;
t91 = (mrSges(5,3) + t85) * t20;
t87 = Ifges(11,3) + Ifges(4,3) + Ifges(9,3) + m(10) * (t38 ^ 2 + t39 ^ 2) * pkin(2) ^ 2;
t45 = -t68 * t69 + t73 * t74;
t34 = -t45 * pkin(5) + t53;
t26 = t38 * t74 + t39 * t69;
t61 = sin(pkin(22));
t64 = cos(pkin(22));
t41 = t61 * t75 - t70 * t64;
t44 = t68 * t74 + t69 * t73;
t60 = qJ(2) + qJ(3);
t55 = sin(t60);
t56 = cos(t60);
t83 = Ifges(11,5) * t55 + Ifges(4,5) * t44 + Ifges(9,5) * t26 + Ifges(11,6) * t56 + Ifges(4,6) * t45 + Ifges(9,6) * t28;
t82 = (mrSges(10,1) * t38 - mrSges(10,2) * t39) * pkin(2);
t13 = t21 * mrSges(6,2) - t67 * t101;
t14 = -t21 * mrSges(6,1) - t72 * t101;
t81 = mrSges(5,3) * t21 + t13 * t72 - t14 * t67;
t76 = cos(pkin(17));
t71 = sin(pkin(17));
t46 = t70 * t71 + t75 * t76;
t43 = t70 * t76 - t71 * t75;
t40 = t61 * t70 + t64 * t75;
t33 = t43 * t69 + t46 * t74;
t32 = t43 * t74 - t46 * t69;
t30 = -t40 * t69 - t41 * t74;
t29 = t40 * t74 - t41 * t69;
t17 = (cos(pkin(21)) * t40 - t41 * sin(pkin(21))) * pkin(4) + t53;
t11 = -t21 * pkin(9) - t20 * pkin(11) + t34;
t6 = t9 ^ 2 * t78;
t1 = t4 ^ 2;
t3 = [Ifges(8,1) * t41 ^ 2 + Ifges(4,2) * t45 ^ 2 + 0.2e1 * (t67 * t13 + t72 * t14) * t11 + (-(mrSges(9,1) * t109) + Ifges(9,2) * t28) * t28 + ((mrSges(9,2) * t109) + Ifges(9,1) * t26 + 0.2e1 * Ifges(9,4) * t28) * t26 + 0.2e1 * (-mrSges(4,1) * t45 + mrSges(8,1) * t40 + mrSges(4,2) * t44 + mrSges(8,2) * t41) * t53 + (-0.2e1 * Ifges(8,4) * t41 + Ifges(8,2) * t40) * t40 + ((mrSges(3,1) * t109) + mrSges(10,1) * t108 + (Ifges(3,2) + Ifges(10,2)) * t69 + 0.2e1 * (-Ifges(3,4) - Ifges(10,4)) * t74) * t69 + ((mrSges(3,2) * t109) + mrSges(10,2) * t108 + (Ifges(10,1) + Ifges(3,1)) * t74) * t74 + (-0.2e1 * mrSges(5,1) * t34 + Ifges(5,2) * t21 - t96) * t21 + m(6) * t95 * t11 ^ 2 + 0.2e1 * t17 * (-mrSges(11,1) * t56 + mrSges(11,2) * t55) + t55 * (Ifges(11,1) * t55 + Ifges(11,4) * t56) + t56 * (Ifges(11,4) * t55 + Ifges(11,2) * t56) + t32 * (Ifges(7,4) * t33 + Ifges(7,2) * t32) + 0.2e1 * pkin(14) * (-mrSges(7,1) * t32 + mrSges(7,2) * t33) + t33 * (Ifges(7,1) * t33 + Ifges(7,4) * t32) + m(5) * t34 ^ 2 + m(10) * t22 ^ 2 + m(11) * t17 ^ 2 + (m(8) + m(4)) * t53 ^ 2 + ((m(9) + m(3)) * pkin(15) ^ 2) + (Ifges(4,1) * t44 + 0.2e1 * Ifges(4,4) * t45) * t44 + m(7) * pkin(14) ^ 2 + (0.2e1 * t34 * mrSges(5,2) + (Ifges(6,1) * t59 + Ifges(5,1) + (-0.2e1 * Ifges(6,4) * t72 + Ifges(6,2) * t67) * t67) * t20 + ((2 * Ifges(5,4)) - t97 + 0.2e1 * t100) * t21) * t20 + Ifges(2,3); Ifges(7,5) * t33 + Ifges(7,6) * t32 + t91 * t2 + (-t38 * t103 + Ifges(3,5) + Ifges(10,5)) * t74 + (-t39 * t103 - Ifges(3,6) - Ifges(10,6)) * t69 + t81 * t4 + ((-t29 * t40 + t30 * t41) * mrSges(8,3) + (-t44 * t68 - t45 * t73) * mrSges(4,3) + (-t55 * t68 - t56 * t73) * mrSges(11,3)) * pkin(1) + t83; Ifges(3,3) + Ifges(7,3) + Ifges(10,3) + m(6) * (t95 * t1 + t110) + m(5) * (t1 + t110) + t87 + (m(8) * (t29 ^ 2 + t30 ^ 2) + 0.2e1 * (m(4) / 0.2e1 + m(11) / 0.2e1) * (t68 ^ 2 + t73 ^ 2)) * pkin(1) ^ 2 + 0.2e1 * t82 + 0.2e1 * t112; -t26 * t103 + (-t111 * t81 - t91 * t9) * pkin(5) + t83; t82 + t112 + 0.2e1 * (m(6) * (-t95 * t105 - t106) / 0.2e1 + m(5) * (-t105 - t106) / 0.2e1) * pkin(5) + t87; m(6) * (t95 * t102 + t6) + m(5) * (t6 + t102) + t87; -t20 * t100 + (mrSges(6,1) * t72 - mrSges(6,2) * t67) * t11 + t96; -t85 * t4; t85 * t111 * pkin(5); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
