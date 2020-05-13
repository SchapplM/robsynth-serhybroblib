% Calculate Gravitation load on the joints for
% fivebar1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fivebar1IC_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1IC_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1IC_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1IC_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:11
% EndTime: 2020-04-27 06:19:11
% DurationCPUTime: 0.12s
% Computational Cost: add. (86->41), mult. (132->54), div. (8->4), fcn. (52->15), ass. (0->28)
t12 = mrSges(5,1) * g(1) + mrSges(5,2) * g(2);
t13 = mrSges(5,1) * g(2) - mrSges(5,2) * g(1);
t18 = sin(qJ(4));
t22 = cos(qJ(4));
t39 = t12 * t18 - t13 * t22;
t14 = mrSges(3,1) * g(1) + mrSges(3,2) * g(2);
t15 = mrSges(3,1) * g(2) - mrSges(3,2) * g(1);
t20 = sin(qJ(2));
t24 = cos(qJ(2));
t38 = -t14 * t24 - t15 * t20;
t19 = sin(qJ(3));
t23 = cos(qJ(3));
t31 = -t12 * t22 - t13 * t18;
t37 = (-t19 * t31 + t39 * t23) / pkin(4);
t21 = sin(qJ(1));
t25 = cos(qJ(1));
t32 = -t14 * t20 + t15 * t24;
t36 = (t21 * t38 + t32 * t25) / pkin(5);
t33 = qJ(4) + qJ(3);
t34 = qJ(1) + qJ(2);
t35 = cos(0.2e1 * t33) - cos(0.2e1 * t34);
t28 = 0.2e1 * qJ(1);
t27 = 0.2e1 * qJ(3);
t17 = pkin(2) * m(3) + mrSges(2,1);
t16 = pkin(3) * m(5) + mrSges(4,1);
t9 = -0.1e1 / sin(t33 - t34);
t3 = 0.1e1 / t35;
t1 = [(mrSges(2,2) * g(1) - g(2) * t17 + t32) * t25 + (mrSges(2,2) * g(2) + t17 * g(1) + t38) * t21 - (t35 * pkin(5) + (-cos(qJ(2) + 0.2e1 * qJ(4) + t27) + cos(qJ(2) + t28)) * pkin(2)) * t3 * t36 + pkin(2) * t20 * t9 * t37; pkin(3) * t18 * t9 * t36 + (mrSges(4,2) * g(1) - g(2) * t16 + t39) * t23 + (mrSges(4,2) * g(2) + t16 * g(1) - t31) * t19 - (t35 * pkin(4) + (cos(qJ(4) + t27) - cos(qJ(4) + t28 + 0.2e1 * qJ(2))) * pkin(3)) * t3 * t37;];
taug = t1(:);
