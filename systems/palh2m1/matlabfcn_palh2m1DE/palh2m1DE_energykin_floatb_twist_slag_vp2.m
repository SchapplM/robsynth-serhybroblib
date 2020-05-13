% Calculate kinetic energy for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m1DE_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:51:54
% EndTime: 2020-05-02 23:51:57
% DurationCPUTime: 2.64s
% Computational Cost: add. (884->232), mult. (1176->308), div. (0->0), fcn. (616->36), ass. (0->102)
t136 = (m(2) * pkin(5));
t57 = V_base(6) + qJD(1);
t135 = t57 / 0.2e1;
t134 = Ifges(2,3) * qJD(1);
t100 = (V_base(5) * V_base(4));
t82 = qJ(1) - qJ(2);
t123 = cos(t82) / 0.2e1;
t124 = sin(t82) / 0.2e1;
t79 = qJ(2) + qJ(3);
t72 = qJ(1) - t79;
t125 = cos(t72) / 0.2e1;
t126 = sin(t72) / 0.2e1;
t71 = qJ(1) + t79;
t44 = sin(t71);
t48 = cos(t71);
t81 = qJ(1) + qJ(2);
t59 = sin(t81);
t63 = cos(t81);
t73 = qJD(2) + qJD(3);
t90 = cos(qJ(2));
t133 = -((t123 + t63 / 0.2e1) * V_base(5) - (t124 + t59 / 0.2e1) * V_base(4) + qJD(2) * t90) * pkin(2) - ((t125 + t48 / 0.2e1) * V_base(5) - (t126 + t44 / 0.2e1) * V_base(4) + t73 * cos(t79)) * pkin(3) + V_base(3);
t127 = -qJD(2) / 0.2e1 + t135;
t87 = sin(qJ(1));
t103 = t87 * V_base(4);
t91 = cos(qJ(1));
t41 = t91 * V_base(5);
t23 = t41 - t103;
t132 = pkin(1) * t23;
t92 = pkin(5) + pkin(6);
t129 = -t92 * V_base(4) + V_base(2);
t128 = -qJD(3) / 0.2e1 + t127;
t85 = sin(qJ(3));
t121 = pkin(3) * t85;
t93 = pkin(1) + pkin(4);
t120 = t57 * t93;
t86 = sin(qJ(2));
t118 = t85 * t86;
t117 = (2 * mrSges(2,3) + t136) * pkin(5);
t75 = V_base(5) ^ 2;
t76 = V_base(4) ^ 2;
t116 = t75 - t76;
t115 = pkin(2) * qJD(2);
t113 = qJ(4) - qJ(2);
t77 = qJ(4) + qJ(2);
t67 = qJ(3) + t77;
t53 = qJ(1) + t67;
t68 = qJ(3) - t113;
t54 = qJ(1) - t68;
t111 = sin(t54) / 0.2e1 - sin(t53) / 0.2e1;
t110 = cos(t54) / 0.2e1 - cos(t53) / 0.2e1;
t69 = qJ(1) + t77;
t70 = qJ(1) + t113;
t109 = sin(t70) / 0.2e1 - sin(t69) / 0.2e1;
t107 = cos(t70) / 0.2e1 - cos(t69) / 0.2e1;
t101 = -V_base(1) * t87 + V_base(2) * t91;
t39 = qJD(2) + t57;
t99 = mrSges(2,1) * t91 - mrSges(2,2) * t87;
t98 = -mrSges(2,1) * t87 - mrSges(2,2) * t91;
t97 = t92 * V_base(5) + V_base(1);
t22 = qJD(2) + t23;
t96 = t87 * V_base(2) + t91 * V_base(1);
t26 = V_base(5) * t87 + V_base(4) * t91;
t89 = cos(qJ(3));
t88 = cos(qJ(4));
t84 = sin(qJ(4));
t83 = -Ifges(2,2) + Ifges(2,1);
t80 = qJ(1) + qJ(4);
t66 = -mrSges(2,1) * pkin(5) + Ifges(2,5);
t65 = mrSges(2,2) * pkin(5) - Ifges(2,6);
t62 = cos(t80);
t58 = sin(t80);
t51 = pkin(3) * t89 + pkin(2);
t38 = qJD(4) + t57;
t37 = -V_base(4) * pkin(5) + V_base(2);
t36 = V_base(5) * pkin(5) + V_base(1);
t30 = qJD(3) + t39;
t29 = pkin(1) * t86 + pkin(5) * t90;
t28 = pkin(1) * t90 - pkin(5) * t86 + pkin(2);
t25 = t89 * t90 - t118;
t24 = t85 * t90 + t86 * t89;
t21 = qJD(3) + t22;
t20 = V_base(3) - t132;
t19 = t36 * t91 + t37 * t87;
t18 = t26 * t90 - t57 * t86;
t17 = -t26 * t86 - t57 * t90;
t16 = pkin(1) * t57 - t36 * t87 + t37 * t91;
t15 = t23 * t88 - t26 * t84;
t14 = t23 * t84 + t26 * t88;
t13 = t19 * t90 - t20 * t86;
t12 = -t19 * t86 - t20 * t90;
t11 = t57 * (pkin(2) * t90 + pkin(1)) + t26 * (pkin(2) * t86 - pkin(5)) + t101;
t10 = t17 * t85 + t18 * t89;
t9 = t17 * t89 - t18 * t85;
t8 = -t24 * V_base(3) + t85 * t115 + t23 * (t28 * t85 + t29 * t89) + t96 * t25;
t7 = -t25 * V_base(3) + t89 * t115 + (t28 * t89 - t29 * t85) * t23 - t96 * t24;
t6 = t132 - t133;
t5 = (-pkin(3) * t118 + t51 * t90 + pkin(1)) * qJD(1) + (t121 * t26 + t51 * V_base(6)) * t90 + (-t121 * V_base(6) + t26 * t51) * t86 + pkin(1) * V_base(6) - t26 * pkin(5) + t101;
t4 = 0.2e1 * (-t41 / 0.2e1 + t103 / 0.2e1) * t93 + t133;
t3 = pkin(5) * t23 + (-t73 * sin(t79) + (t126 - t44 / 0.2e1) * V_base(5) + (t125 - t48 / 0.2e1) * V_base(4)) * pkin(3) + (-qJD(2) * t86 + (t124 - t59 / 0.2e1) * V_base(5) + (t123 - t63 / 0.2e1) * V_base(4)) * pkin(2) + t96;
t2 = t97 * t62 + t58 * t129 + t84 * t120 + (-t30 * sin(t68) / 0.2e1 + sin(t67) * t128 + t111 * V_base(5) + t110 * V_base(4)) * pkin(3) + (t39 * sin(t113) / 0.2e1 + sin(t77) * t127 + t109 * V_base(5) + t107 * V_base(4)) * pkin(2);
t1 = t62 * t129 - t97 * t58 + t88 * t120 + (t30 * cos(t68) / 0.2e1 + cos(t67) * t128 + t110 * V_base(5) - t111 * V_base(4)) * pkin(3) + (t39 * cos(t113) / 0.2e1 + cos(t77) * t127 + t107 * V_base(5) - t109 * V_base(4)) * pkin(2);
t27 = ((Ifges(1,4) - Ifges(2,4)) * t100) + ((4 * Ifges(2,4) * t100) - t116 * t83) * t91 ^ 2 / 0.2e1 + (Ifges(2,1) + Ifges(1,2) + t117) * t75 / 0.2e1 + (Ifges(1,1) + Ifges(2,2) + t117) * t76 / 0.2e1 + m(4) * (t11 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + m(5) * (t3 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (-V_base(5) * mrSges(1,1) + V_base(4) * mrSges(1,2) + (-V_base(5) * mrSges(2,1) + V_base(4) * mrSges(2,2)) * t91 + (V_base(4) * mrSges(2,1) + V_base(5) * mrSges(2,2)) * t87) * V_base(3) + (t4 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,4) * t15 + Ifges(6,5) * t38 + Ifges(6,1) * t14 / 0.2e1) * t14 + (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) * (m(1) / 0.2e1 + m(2) / 0.2e1) + (t12 * mrSges(3,1) - t13 * mrSges(3,2) + Ifges(3,3) * t22 / 0.2e1) * t22 + ((-t65 * t91 + t66 * t87) * V_base(5) + t98 * V_base(1) + t99 * V_base(2) + (t65 * t87 + t66 * t91) * V_base(4) + t134 / 0.2e1) * qJD(1) + (V_base(5) * Ifges(1,6) + (t65 * V_base(4) + t66 * V_base(5)) * t87 + Ifges(1,5) * V_base(4) + (-t65 * V_base(5) + t66 * V_base(4)) * t91 + t134 + (-mrSges(1,2) + t98) * V_base(1) + (mrSges(1,1) + t99) * V_base(2) + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (t5 * mrSges(5,2) - t6 * mrSges(5,3) - Ifges(5,4) * t57 + Ifges(5,1) * t26 / 0.2e1) * t26 + (-t4 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,6) * t38 + Ifges(6,2) * t15 / 0.2e1) * t15 + (t7 * mrSges(4,1) - t8 * mrSges(4,2) + Ifges(4,6) * t9 + Ifges(4,3) * t21 / 0.2e1) * t21 + (-t11 * mrSges(4,1) + t8 * mrSges(4,3) + Ifges(4,2) * t9 / 0.2e1) * t9 + (t11 * mrSges(4,2) - t7 * mrSges(4,3) + Ifges(4,4) * t9 + Ifges(4,5) * t21 + Ifges(4,1) * t10 / 0.2e1) * t10 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t38 / 0.2e1) * t38 + (t6 * mrSges(5,1) - t3 * mrSges(5,2) + Ifges(5,5) * t26 - Ifges(5,6) * t57 + Ifges(5,3) * t23 / 0.2e1) * t23 + (-t16 * mrSges(3,1) + t13 * mrSges(3,3) + Ifges(3,4) * t18 + Ifges(3,6) * t22 + Ifges(3,2) * t17 / 0.2e1) * t17 + m(3) * (t12 ^ 2 + t13 ^ 2 + t16 ^ 2) / 0.2e1 + (Ifges(2,4) * t116 + t100 * t83) * t87 * t91 + (t5 * mrSges(5,1) - t3 * mrSges(5,3) + Ifges(5,2) * t135) * t57 + (t16 * mrSges(3,2) - t12 * mrSges(3,3) + Ifges(3,5) * t22 + Ifges(3,1) * t18 / 0.2e1) * t18 + (-V_base(4) * V_base(2) + V_base(5) * V_base(1)) * (mrSges(1,3) + mrSges(2,3) + t136);
T = t27;
