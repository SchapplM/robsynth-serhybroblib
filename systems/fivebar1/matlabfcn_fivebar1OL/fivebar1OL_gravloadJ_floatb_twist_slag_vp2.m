% Calculate Gravitation load on the joints for
% fivebar1OL
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
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fivebar1OL_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1OL_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:00
% EndTime: 2020-04-27 06:13:01
% DurationCPUTime: 0.09s
% Computational Cost: add. (40->22), mult. (68->34), div. (0->0), fcn. (24->8), ass. (0->19)
t12 = sin(qJ(4));
t16 = cos(qJ(4));
t6 = mrSges(5,1) * g(1) + mrSges(5,2) * g(2);
t7 = mrSges(5,1) * g(2) - mrSges(5,2) * g(1);
t24 = t6 * t12 - t7 * t16;
t14 = sin(qJ(2));
t18 = cos(qJ(2));
t8 = mrSges(3,1) * g(1) + mrSges(3,2) * g(2);
t9 = mrSges(3,1) * g(2) - mrSges(3,2) * g(1);
t23 = -t9 * t14 - t8 * t18;
t22 = -t8 * t14 + t9 * t18;
t21 = -t7 * t12 - t6 * t16;
t19 = cos(qJ(1));
t17 = cos(qJ(3));
t15 = sin(qJ(1));
t13 = sin(qJ(3));
t11 = pkin(2) * m(3) + mrSges(2,1);
t10 = pkin(3) * m(5) + mrSges(4,1);
t1 = [(mrSges(2,2) * g(1) - g(2) * t11 + t22) * t19 + (mrSges(2,2) * g(2) + t11 * g(1) + t23) * t15, t15 * t23 + t22 * t19, (mrSges(4,2) * g(1) - g(2) * t10 + t24) * t17 + (mrSges(4,2) * g(2) + t10 * g(1) - t21) * t13, -t13 * t21 + t24 * t17, 0];
taug = t1(:);
