% Calculate joint inertia matrix for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m2DE_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2DE_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2DE_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:26
% EndTime: 2020-05-03 01:06:26
% DurationCPUTime: 0.10s
% Computational Cost: add. (88->46), mult. (107->64), div. (0->0), fcn. (62->6), ass. (0->19)
t15 = m(6) + m(7);
t11 = sin(qJ(2));
t20 = pkin(4) * t11;
t10 = sin(qJ(3));
t19 = pkin(5) * t10;
t18 = t10 * mrSges(5,2);
t14 = cos(qJ(2));
t7 = -t14 * pkin(4) - pkin(1);
t12 = cos(qJ(4));
t9 = sin(qJ(4));
t5 = mrSges(7,1) * t9 + mrSges(7,2) * t12;
t6 = -pkin(2) + t7;
t17 = -mrSges(6,3) + t5;
t13 = cos(qJ(3));
t4 = -t13 * pkin(5) + t6;
t2 = pkin(3) - t4;
t16 = (mrSges(7,1) * t12 - mrSges(7,2) * t9) * t2;
t8 = t15 * pkin(5) + mrSges(5,1);
t1 = [Ifges(4,2) + Ifges(6,2) + Ifges(2,3) + Ifges(7,3) + m(5) * t6 ^ 2 + m(3) * pkin(1) ^ 2 + 0.2e1 * t6 * (-t13 * mrSges(5,1) + t18) + Ifges(5,2) * t13 ^ 2 - 0.2e1 * pkin(1) * (-t14 * mrSges(3,1) + t11 * mrSges(3,2)) + t11 * (Ifges(3,1) * t11 + Ifges(3,4) * t14) + t14 * (Ifges(3,4) * t11 + Ifges(3,2) * t14) + (m(4) * t7 - (2 * mrSges(4,1))) * t7 + (m(6) * t4 - (2 * mrSges(6,1))) * t4 + 0.2e1 * t16 + m(7) * (t12 ^ 2 + t9 ^ 2) * t2 ^ 2 + (Ifges(5,1) * t10 + 0.2e1 * Ifges(5,4) * t13) * t10; Ifges(3,5) * t11 + Ifges(3,6) * t14 + (-mrSges(4,3) - mrSges(5,3) + t17) * t20; (m(4) + m(5) + t15) * pkin(4) ^ 2 + Ifges(3,3); Ifges(5,5) * t10 + Ifges(5,6) * t13 + t17 * t19; pkin(4) * ((t8 * t13 - t18) * t14 + (mrSges(5,2) * t13 + t8 * t10) * t11); t15 * pkin(5) ^ 2 + Ifges(5,3); Ifges(7,3) + t16; t5 * t20; t5 * t19; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
