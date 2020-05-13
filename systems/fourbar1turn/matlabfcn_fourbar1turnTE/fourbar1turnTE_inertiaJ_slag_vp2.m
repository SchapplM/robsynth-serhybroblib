% Calculate joint inertia matrix for
% fourbar1turnTE
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
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnTE_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_inertiaJ_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnTE_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnTE_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:36
% EndTime: 2020-04-12 19:18:38
% DurationCPUTime: 0.52s
% Computational Cost: add. (3841->72), mult. (5259->153), div. (232->12), fcn. (1466->4), ass. (0->56)
t73 = pkin(3) ^ 2;
t30 = sin(qJ(2));
t29 = t30 ^ 2;
t37 = pkin(2) ^ 2;
t38 = pkin(1) ^ 2;
t31 = cos(qJ(2));
t61 = pkin(2) * t31;
t49 = -0.2e1 * t61;
t51 = pkin(1) * t49 + t38;
t25 = t37 + t51;
t50 = -pkin(4) ^ 2 + t73;
t22 = t25 - t50;
t26 = pkin(1) - t61;
t67 = -pkin(3) - pkin(4);
t19 = (pkin(2) - t67) * (pkin(2) + t67) + t51;
t66 = -pkin(3) + pkin(4);
t20 = (pkin(2) - t66) * (pkin(2) + t66) + t51;
t39 = sqrt(-t19 * t20);
t53 = t30 * t39;
t11 = -pkin(2) * t53 + t26 * t22;
t70 = -t11 / 0.2e1;
t62 = pkin(2) * t30;
t17 = t22 * t62;
t12 = t26 * t39 + t17;
t69 = t12 / 0.2e1;
t72 = -0.2e1 * pkin(2);
t68 = -t31 / 0.2e1;
t65 = m(4) * t37;
t60 = t30 * pkin(1);
t48 = pkin(2) * t60;
t57 = 0.1e1 / t39 * (-t19 - t20) * t48;
t23 = 0.1e1 / t25;
t33 = 0.1e1 / pkin(4);
t56 = t23 * t33;
t35 = 0.1e1 / pkin(3);
t55 = t23 * t35;
t52 = t31 * t39;
t47 = pkin(2) * t55;
t24 = 0.1e1 / t25 ^ 2;
t46 = t24 * t62;
t21 = t25 + t50;
t27 = pkin(1) * t31 - pkin(2);
t10 = -pkin(1) * t53 - t27 * t21;
t42 = t10 * t47;
t18 = t21 * t60;
t13 = -t27 * t39 + t18;
t41 = t13 * t47;
t40 = t10 ^ 2;
t9 = t13 ^ 2;
t8 = 0.1e1 / t11 ^ 2;
t7 = 0.1e1 / t40;
t4 = (t30 * t13 / 0.2e1 + t10 * t68) * t55;
t3 = (-t30 * t10 / 0.2e1 + t13 * t68) * t55;
t2 = 0.1e1 + (((0.2e1 * t38 * t29 * pkin(2) - t27 * t57) * t23 + ((t31 * t21 + t53) * t23 - 0.2e1 * t13 * t46) * pkin(1)) / t10 - (t18 * t23 + (-t23 * t52 + ((t27 * t72 - t57) * t23 + t10 * t24 * t72) * t30) * pkin(1)) * t13 * t7) * pkin(3) * t25 * t35 / (t9 * t7 + 0.1e1);
t1 = 0.2e1 * (-((t26 * t57 + (t31 * t22 + t53) * pkin(2)) * t23 / 0.2e1 + (t37 * t29 * t23 - t12 * t46) * pkin(1)) / t11 - (-(t17 + (-t30 * t57 - t52) * pkin(2)) * t23 / 0.2e1 + (t11 * t24 - t23 * t26) * t48) * t12 * t8) * pkin(4) * t25 * t33 / (t12 ^ 2 * t8 + 0.1e1);
t5 = [m(5) * t38 - 0.2e1 * pkin(1) * (mrSges(5,2) * t69 + t11 * mrSges(5,1) / 0.2e1) * t56 + Ifges(3,1) * t29 + Ifges(2,3) + (-mrSges(4,1) * t49 + Ifges(4,2) * t4) * t4 + (mrSges(4,2) * t49 + Ifges(4,1) * t3 + 0.2e1 * Ifges(4,4) * t4) * t3 + (Ifges(5,2) * t70 ^ 2 + (Ifges(5,1) * t69 + 0.2e1 * Ifges(5,4) * t70) * t69) * t23 ^ 2 * t33 ^ 2 + (0.2e1 * Ifges(3,4) * t30 + (t65 + Ifges(3,2)) * t31) * t31; Ifges(3,5) * t30 + Ifges(3,6) * t31 + t2 * (Ifges(4,5) * t3 + Ifges(4,6) * t4) + (-t4 * t41 / 0.2e1 + t3 * t42 / 0.2e1) * mrSges(4,3) + (Ifges(5,5) * t69 + Ifges(5,6) * t70) * t1 * t56; t1 ^ 2 * Ifges(5,3) + Ifges(3,3) + (t9 / 0.2e1 + t40 / 0.2e1) / t73 * t24 * t65 / 0.2e1 + (-mrSges(4,1) * t42 + mrSges(4,2) * t41 + Ifges(4,3) * t2) * t2;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t5(1), t5(2); t5(2), t5(3);];
Mq = res;
