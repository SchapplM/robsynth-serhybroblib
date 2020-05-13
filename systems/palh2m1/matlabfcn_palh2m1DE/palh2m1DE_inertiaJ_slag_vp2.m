% Calculate joint inertia matrix for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m1DE_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1DE_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1DE_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:04
% EndTime: 2020-05-02 23:52:05
% DurationCPUTime: 0.22s
% Computational Cost: add. (195->87), mult. (227->89), div. (0->0), fcn. (72->22), ass. (0->43)
t25 = qJ(4) + qJ(2);
t20 = qJ(3) + t25;
t7 = sin(t20);
t26 = -qJ(4) + qJ(2);
t21 = qJ(3) + t26;
t8 = sin(t21);
t58 = (t7 + t8) * mrSges(6,2);
t10 = cos(t21);
t9 = cos(t20);
t51 = -t9 / 0.2e1;
t57 = (t10 / 0.2e1 + t51) * mrSges(6,1);
t34 = cos(qJ(2));
t56 = 0.2e1 * t34;
t36 = m(5) + m(6);
t11 = t36 * pkin(3) + mrSges(4,1);
t14 = sin(t25);
t15 = sin(t26);
t17 = cos(t25);
t18 = cos(t26);
t55 = (t15 / 0.2e1 + t14 / 0.2e1) * mrSges(6,2) + (t18 / 0.2e1 - t17 / 0.2e1) * mrSges(6,1);
t22 = (m(6) * pkin(4) + mrSges(5,1));
t53 = 2 * t22;
t52 = pkin(3) / 0.2e1;
t49 = pkin(1) + pkin(4);
t38 = pkin(3) ^ 2;
t23 = t36 * t38;
t46 = Ifges(4,3) + t23;
t24 = m(4) + t36;
t30 = sin(qJ(3));
t45 = t30 * pkin(2) * mrSges(4,2);
t41 = mrSges(6,1) * cos(qJ(4)) - mrSges(6,2) * sin(qJ(4));
t28 = qJ(2) + qJ(3);
t16 = sin(t28);
t19 = cos(t28);
t40 = -Ifges(4,6) * t19 + (mrSges(5,3) * pkin(3) - Ifges(4,5)) * t16 + pkin(3) * t57 + t52 * t58;
t37 = 0.2e1 * qJ(2);
t33 = cos(qJ(3));
t31 = sin(qJ(2));
t27 = t37 + qJ(3);
t13 = t24 * pkin(2) ^ 2;
t12 = 0.2e1 * t28;
t2 = t11 * pkin(2) * t33;
t1 = [(t23 - Ifges(4,1) + Ifges(4,2)) * cos(t12) / 0.2e1 + Ifges(3,4) * sin(t37) + (t13 - Ifges(3,1) + Ifges(3,2)) * cos(t37) / 0.2e1 + Ifges(4,4) * sin(t12) + t13 / 0.2e1 + m(6) * pkin(4) ^ 2 + t2 + Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,2) + Ifges(2,3) + Ifges(6,3) + t23 / 0.2e1 + 0.2e1 * t41 * t49 + (t19 * t53 + (-t7 + t8) * mrSges(6,2) + (t10 + t9) * mrSges(6,1)) * pkin(3) + (mrSges(3,1) * t56 - 0.2e1 * mrSges(4,2) * t16 - 0.2e1 * mrSges(3,2) * t31 + (m(3) + t24) * pkin(1) + t53 + 0.2e1 * t11 * t19) * pkin(1) + ((pkin(1) * t24 + t22) * t56 + t11 * cos(t27) + (-t14 + t15) * mrSges(6,2) + (-sin(t27) - t30) * mrSges(4,2) + (t17 + t18) * mrSges(6,1)) * pkin(2); -Ifges(3,5) * t31 - Ifges(3,6) * t34 + ((mrSges(4,3) + mrSges(5,3)) * t31 + t55) * pkin(2) + t40; Ifges(3,3) + t13 + 0.2e1 * t2 - 0.2e1 * t45 + t46; t40; t2 - t45 + t46; Ifges(4,3) + (m(5) + (t10 * t51 + t8 * t7 / 0.2e1 + t19 ^ 2 + 0.1e1 / 0.2e1) * m(6)) * t38; Ifges(6,3) + t41 * ((t33 * pkin(3) + pkin(2)) * t34 - t31 * t30 * pkin(3) + t49); (-(-t8 / 0.2e1 - t7 / 0.2e1) * mrSges(6,2) + t57) * pkin(3) + t55 * pkin(2); (t58 + (t10 - t9) * mrSges(6,1)) * t52; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
