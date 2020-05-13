% Calculate kinetic energy for
% palh1m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m1OL_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(6,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:17
% EndTime: 2020-04-15 19:28:23
% DurationCPUTime: 4.42s
% Computational Cost: add. (2913->238), mult. (3985->358), div. (0->0), fcn. (3216->22), ass. (0->98)
t107 = cos(qJ(1));
t98 = sin(qJ(1));
t77 = t107 * V_base(5) - t98 * V_base(4);
t74 = qJD(2) - t77;
t71 = qJD(7) + t74;
t110 = pkin(4) * t71;
t104 = cos(qJ(4));
t105 = cos(qJ(3));
t106 = cos(qJ(2));
t82 = V_base(5) * pkin(13) + V_base(1);
t83 = -V_base(4) * pkin(13) + V_base(2);
t64 = t107 * t82 + t98 * t83;
t66 = -pkin(15) * t77 + V_base(3);
t97 = sin(qJ(2));
t54 = -t106 * t64 - t66 * t97;
t39 = pkin(1) * t74 + t54;
t53 = t106 * t66 - t64 * t97;
t96 = sin(qJ(3));
t31 = t105 * t53 + t96 * t39;
t72 = qJD(3) + t74;
t20 = pkin(5) * t72 + t31;
t29 = -t105 * t39 + t53 * t96;
t95 = sin(qJ(4));
t9 = t104 * t29 + t95 * t20;
t101 = cos(qJ(7));
t92 = sin(qJ(7));
t30 = t101 * t53 + t92 * t39;
t109 = sin(qJ(10));
t28 = t101 * t39 - t53 * t92;
t100 = cos(qJ(8));
t91 = sin(qJ(8));
t33 = t100 * t54 - t53 * t91;
t63 = t107 * t83 - t98 * t82;
t8 = t104 * t20 - t29 * t95;
t78 = t107 * V_base(4) + t98 * V_base(5);
t86 = V_base(6) + qJD(1);
t60 = t106 * t86 - t78 * t97;
t62 = -t106 * t78 - t86 * t97;
t44 = -t105 * t62 + t60 * t96;
t47 = t105 * t60 + t62 * t96;
t23 = t104 * t47 - t44 * t95;
t58 = -pkin(15) * t86 - t63;
t70 = qJD(8) + t74;
t49 = -pkin(1) * t62 + t58;
t32 = -pkin(5) * t47 + t49;
t108 = V_base(3) ^ 2;
t103 = cos(qJ(5));
t102 = cos(qJ(6));
t99 = cos(qJ(9));
t94 = sin(qJ(5));
t93 = sin(qJ(6));
t90 = sin(qJ(9));
t89 = cos(qJ(10));
t88 = cos(pkin(19));
t87 = sin(pkin(19));
t76 = -t87 * t109 - t88 * t89;
t75 = -t88 * t109 + t87 * t89;
t73 = qJD(6) - t77;
t69 = qJD(4) + t72;
t68 = qJD(9) + t70;
t67 = pkin(14) * t77 + V_base(3);
t65 = qJD(10) + t71;
t61 = t102 * t78 + t86 * t93;
t59 = t102 * t86 - t78 * t93;
t57 = t58 ^ 2;
t56 = -pkin(16) * t77 + t64;
t55 = pkin(14) * t86 - pkin(16) * t78 - t63;
t48 = t49 ^ 2;
t46 = t101 * t60 + t62 * t92;
t45 = t100 * t60 + t62 * t91;
t43 = t101 * t62 - t60 * t92;
t42 = t100 * t62 - t60 * t91;
t41 = t102 * t56 + t67 * t93;
t40 = t102 * t67 - t56 * t93;
t35 = -pkin(2) * t42 + t58;
t34 = t100 * t53 + t54 * t91;
t27 = pkin(2) * t70 + t33;
t25 = t104 * t44 + t47 * t95;
t24 = -t42 * t90 - t45 * t99;
t22 = -t42 * t99 + t45 * t90;
t21 = qJD(5) - t23;
t19 = t87 * t110 + t30;
t18 = t88 * t110 + t28;
t16 = t43 * t75 + t46 * t76;
t15 = t43 * t76 - t46 * t75;
t14 = t103 * t25 + t69 * t94;
t13 = t103 * t69 - t25 * t94;
t12 = (-t43 * t88 - t46 * t87) * pkin(4) + t49;
t11 = -t27 * t90 - t34 * t99;
t10 = -t27 * t99 + t34 * t90;
t7 = pkin(11) * t69 + t9;
t6 = -pkin(9) * t69 - t8;
t5 = t18 * t75 + t19 * t76;
t4 = t18 * t76 - t19 * t75;
t3 = -pkin(9) * t23 - pkin(11) * t25 + t32;
t2 = t103 * t7 + t3 * t94;
t1 = t103 * t3 - t7 * t94;
t17 = (-t12 * mrSges(11,1) + t5 * mrSges(11,3) + Ifges(11,4) * t16 + Ifges(11,6) * t65 + Ifges(11,2) * t15 / 0.2e1) * t15 + (-t6 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t14 + Ifges(6,6) * t21 + Ifges(6,2) * t13 / 0.2e1) * t13 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t55 * mrSges(7,1) + t41 * mrSges(7,3) + Ifges(7,4) * t61 + Ifges(7,6) * t73 + Ifges(7,2) * t59 / 0.2e1) * t59 + (t35 * mrSges(10,2) - t10 * mrSges(10,3) + Ifges(10,5) * t68 + Ifges(10,1) * t24 / 0.2e1) * t24 + (-t58 * mrSges(3,1) + t53 * mrSges(3,3) + Ifges(3,6) * t74 + Ifges(3,2) * t62 / 0.2e1) * t62 + (-V_base(3) * mrSges(2,1) + t64 * mrSges(2,3) + Ifges(2,4) * t78 + Ifges(2,6) * t86 + Ifges(2,2) * t77 / 0.2e1) * t77 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-t49 * mrSges(4,1) + t29 * mrSges(4,3) + Ifges(4,6) * t72 + Ifges(4,2) * t47 / 0.2e1) * t47 + (t4 * mrSges(11,1) - t5 * mrSges(11,2) + Ifges(11,3) * t65 / 0.2e1) * t65 + (t10 * mrSges(10,1) - t11 * mrSges(10,2) + Ifges(10,3) * t68 / 0.2e1) * t68 + (t49 * mrSges(8,2) - t28 * mrSges(8,3) + Ifges(8,5) * t71 + Ifges(8,1) * t46 / 0.2e1) * t46 + m(2) * (t63 ^ 2 + t64 ^ 2 + t108) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t108) / 0.2e1 + (t54 * mrSges(3,1) - t53 * mrSges(3,2) + Ifges(3,3) * t74 / 0.2e1) * t74 + (t40 * mrSges(7,1) - t41 * mrSges(7,2) + Ifges(7,3) * t73 / 0.2e1) * t73 + (t49 * mrSges(4,2) - t31 * mrSges(4,3) + Ifges(4,4) * t47 + Ifges(4,5) * t72 + Ifges(4,1) * t44 / 0.2e1) * t44 + (t58 * mrSges(9,2) - t33 * mrSges(9,3) + Ifges(9,5) * t70 + Ifges(9,1) * t45 / 0.2e1) * t45 + (t33 * mrSges(9,1) - t34 * mrSges(9,2) + Ifges(9,3) * t70 / 0.2e1) * t70 + (-t49 * mrSges(8,1) + t30 * mrSges(8,3) + Ifges(8,4) * t46 + Ifges(8,6) * t71 + Ifges(8,2) * t43 / 0.2e1) * t43 + (-t58 * mrSges(9,1) + t34 * mrSges(9,3) + Ifges(9,4) * t45 + Ifges(9,6) * t70 + Ifges(9,2) * t42 / 0.2e1) * t42 + (t28 * mrSges(8,1) - t30 * mrSges(8,2) + Ifges(8,3) * t71 / 0.2e1) * t71 + (t6 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t21 + Ifges(6,1) * t14 / 0.2e1) * t14 + m(4) * (t29 ^ 2 + t31 ^ 2 + t48) / 0.2e1 + m(7) * (t40 ^ 2 + t41 ^ 2 + t55 ^ 2) / 0.2e1 + m(3) * (t53 ^ 2 + t54 ^ 2 + t57) / 0.2e1 + m(9) * (t33 ^ 2 + t34 ^ 2 + t57) / 0.2e1 + m(8) * (t28 ^ 2 + t30 ^ 2 + t48) / 0.2e1 + m(10) * (t10 ^ 2 + t11 ^ 2 + t35 ^ 2) / 0.2e1 + m(5) * (t32 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(11) * (t12 ^ 2 + t4 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2) / 0.2e1 + (t32 * mrSges(5,2) - t8 * mrSges(5,3) + Ifges(5,5) * t69 + Ifges(5,1) * t25 / 0.2e1) * t25 + (t31 * mrSges(4,1) - t29 * mrSges(4,2) + Ifges(4,3) * t72 / 0.2e1) * t72 + (-t32 * mrSges(5,1) + t9 * mrSges(5,3) + Ifges(5,4) * t25 + Ifges(5,6) * t69 + Ifges(5,2) * t23 / 0.2e1) * t23 + (t55 * mrSges(7,2) - t40 * mrSges(7,3) + Ifges(7,5) * t73 + Ifges(7,1) * t61 / 0.2e1) * t61 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t58 * mrSges(3,2) - t54 * mrSges(3,3) + Ifges(3,4) * t62 + Ifges(3,5) * t74 + Ifges(3,1) * t60 / 0.2e1) * t60 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t21 / 0.2e1) * t21 + (V_base(3) * mrSges(2,2) - t63 * mrSges(2,3) + Ifges(2,5) * t86 + Ifges(2,1) * t78 / 0.2e1) * t78 + (-t35 * mrSges(10,1) + t11 * mrSges(10,3) + Ifges(10,4) * t24 + Ifges(10,6) * t68 + Ifges(10,2) * t22 / 0.2e1) * t22 + (t63 * mrSges(2,1) - t64 * mrSges(2,2) + Ifges(2,3) * t86 / 0.2e1) * t86 + (t12 * mrSges(11,2) - t4 * mrSges(11,3) + Ifges(11,5) * t65 + Ifges(11,1) * t16 / 0.2e1) * t16 + (t8 * mrSges(5,1) - t9 * mrSges(5,2) + Ifges(5,3) * t69 / 0.2e1) * t69;
T = t17;
