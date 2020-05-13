% Calculate Gravitation load on the joints for
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% m [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [1x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1DE1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1DE1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1DE1_gravloadJ_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1DE1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:57:10
% EndTime: 2020-04-24 19:57:12
% DurationCPUTime: 0.30s
% Computational Cost: add. (325->53), mult. (502->81), div. (16->5), fcn. (140->4), ass. (0->35)
t25 = cos(qJ(1));
t41 = pkin(2) * t25;
t37 = (-0.2e1 * t41 + pkin(1)) * pkin(1);
t44 = -pkin(3) - pkin(4);
t4 = (pkin(2) - t44) * (pkin(2) + t44) + t37;
t43 = -pkin(3) + pkin(4);
t5 = (pkin(2) - t43) * (pkin(2) + t43) + t37;
t32 = sqrt(-t4 * t5);
t24 = sin(qJ(1));
t40 = t24 * pkin(2);
t35 = pkin(1) * t40;
t46 = (-t4 - t5) * t35 / t32;
t18 = pkin(2) ^ 2 + t37;
t16 = 0.1e1 / t18;
t45 = t16 / 0.2e1;
t29 = 0.1e1 / pkin(3);
t42 = pkin(1) * t29;
t39 = mrSges(3,2) * t24;
t38 = mrSges(4,2) * t24;
t12 = mrSges(4,1) * t40 + mrSges(4,2) * t41;
t14 = mrSges(3,1) * t40 + mrSges(3,2) * t41;
t36 = pkin(3) ^ 2 - pkin(4) ^ 2;
t34 = t29 * t45;
t33 = pkin(1) - t41;
t27 = 0.1e1 / pkin(4);
t17 = 0.1e1 / t18 ^ 2;
t15 = (mrSges(3,1) * t25 - t39) * pkin(2);
t13 = (mrSges(4,1) * t25 - t38) * pkin(2);
t11 = t18 - t36;
t10 = t18 + t36;
t9 = mrSges(3,1) * t33 + pkin(2) * t39;
t8 = -mrSges(3,2) * pkin(1) + t14;
t7 = -mrSges(4,2) * pkin(1) + t12;
t6 = mrSges(4,1) * t33 + pkin(2) * t38;
t1 = [(mrSges(2,2) * t24 - mrSges(2,1) * t25 - (-t10 * t15 + t14 * t32 + t46 * t9) * t34 - (m(3) * t25 + (-t8 * t16 - (-t10 * t8 + t32 * t9) * t17) * t24 * t42) * pkin(2) - ((t11 * t13 + t12 * t32 + t46 * t6) * t45 + (t7 * t16 - (t11 * t7 + t32 * t6) * t17) * t35) * t27) * g(2) + (mrSges(2,1) * t24 + mrSges(2,2) * t25 - (t10 * t14 + t15 * t32 + t46 * t8) * t34 - (-m(3) + (t9 * t16 - (t10 * t9 + t32 * t8) * t17) * t42) * t40 - ((-t11 * t12 + t13 * t32 + t46 * t7) * t45 + (-t6 * t16 - (-t11 * t6 + t32 * t7) * t17) * t35) * t27) * g(1)];
taug = t1(:);
