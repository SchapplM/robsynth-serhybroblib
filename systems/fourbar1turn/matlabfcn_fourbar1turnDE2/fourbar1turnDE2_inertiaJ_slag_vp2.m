% Calculate joint inertia matrix for
% fourbar1turnDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnDE2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_inertiaJ_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE2_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnDE2_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:31
% EndTime: 2020-04-12 19:33:34
% DurationCPUTime: 0.60s
% Computational Cost: add. (5819->73), mult. (8161->143), div. (418->13), fcn. (2290->8), ass. (0->53)
t44 = pkin(2) ^ 2;
t45 = pkin(1) ^ 2;
t37 = cos(qJ(2));
t60 = t37 * pkin(1);
t66 = -2 * pkin(2);
t56 = t60 * t66 + t45;
t31 = t44 + t56;
t67 = pkin(3) ^ 2;
t68 = pkin(4) ^ 2;
t55 = t67 - t68;
t27 = t31 + t55;
t33 = -pkin(2) + t60;
t36 = sin(qJ(2));
t65 = (-pkin(3) - pkin(4));
t25 = ((pkin(2) - t65) * (pkin(2) + t65)) + t56;
t64 = (-pkin(3) + pkin(4));
t26 = ((pkin(2) - t64) * (pkin(2) + t64)) + t56;
t46 = sqrt(-t25 * t26);
t58 = t36 * t46;
t16 = -pkin(1) * t58 - t27 * t33;
t69 = t16 ^ 2;
t28 = t31 - t55;
t62 = pkin(2) * t36;
t23 = t28 * t62;
t32 = -pkin(2) * t37 + pkin(1);
t18 = t32 * t46 + t23;
t14 = t18 ^ 2;
t35 = t36 ^ 2;
t63 = m(4) * t44;
t61 = t36 * pkin(1);
t54 = pkin(2) * t61;
t59 = 0.1e1 / t46 * (-t25 - t26) * t54;
t57 = t37 * t46;
t24 = t27 * t61;
t19 = -t33 * t46 + t24;
t15 = t19 ^ 2;
t29 = 0.1e1 / t31;
t30 = 0.1e1 / t31 ^ 2;
t42 = 0.1e1 / pkin(3);
t53 = t29 * t42 * ((t15 + t69) / t67 * t30) ^ (-0.1e1 / 0.2e1);
t52 = t30 * t62;
t39 = 0.1e1 / pkin(4);
t17 = -pkin(2) * t58 + t28 * t32;
t47 = t17 ^ 2;
t51 = t39 * t29 * ((t14 + t47) / t68 * t30) ^ (-0.1e1 / 0.2e1);
t50 = pkin(2) * t53;
t13 = 0.1e1 / t47;
t12 = 0.1e1 / t69;
t4 = (-t16 * t37 + t19 * t36) * t53;
t3 = (-t16 * t36 - t19 * t37) * t53;
t2 = 0.1e1 + (((0.2e1 * t45 * t35 * pkin(2) - t33 * t59) * t29 + ((t37 * t27 + t58) * t29 - 0.2e1 * t19 * t52) * pkin(1)) / t16 - (t24 * t29 + (-t29 * t57 + ((t33 * t66 - t59) * t29 + t16 * t30 * t66) * t36) * pkin(1)) * t19 * t12) * pkin(3) / (t12 * t15 + 0.1e1) * t31 * t42;
t1 = 0.2e1 * (-((t32 * t59 + (t37 * t28 + t58) * pkin(2)) * t29 / 0.2e1 + (t29 * t35 * t44 - t18 * t52) * pkin(1)) / t17 - (-(t23 + (-t36 * t59 - t57) * pkin(2)) * t29 / 0.2e1 + (t17 * t30 - t29 * t32) * t54) * t18 * t13) * pkin(4) * t31 * t39 / (t13 * t14 + 0.1e1);
t5 = [Ifges(3,1) * t35 + Ifges(4,2) * t4 ^ 2 + m(5) * t45 + Ifges(2,3) + (Ifges(4,1) * t3 + 0.2e1 * Ifges(4,4) * t4) * t3 + ((-mrSges(4,1) * t4 + mrSges(4,2) * t3) * t66 + 0.2e1 * Ifges(3,4) * t36 + (Ifges(3,2) + t63) * t37) * t37 + (-0.2e1 * pkin(1) * (mrSges(5,1) * t17 + mrSges(5,2) * t18) + (Ifges(5,1) * t14 + (-0.2e1 * Ifges(5,4) * t18 + Ifges(5,2) * t17) * t17) * t51) * t51; Ifges(3,5) * t36 + Ifges(3,6) * t37 + (t16 * t3 - t19 * t4) * mrSges(4,3) * t50 + t2 * (Ifges(4,5) * t3 + Ifges(4,6) * t4) + t1 * (Ifges(5,5) * t18 - Ifges(5,6) * t17) * t51; t1 ^ 2 * Ifges(5,3) + Ifges(3,3) + t63 + (Ifges(4,3) * t2 + 0.2e1 * (-mrSges(4,1) * t16 + mrSges(4,2) * t19) * t50) * t2;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t5(1), t5(2); t5(2), t5(3);];
Mq = res;
