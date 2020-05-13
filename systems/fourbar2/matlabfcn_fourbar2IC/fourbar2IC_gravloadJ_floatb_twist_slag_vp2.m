% Calculate Gravitation load on the joints for
% fourbar2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
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
% Datum: 2020-04-24 20:37
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar2IC_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(2,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2IC_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2IC_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2IC_gravloadJ_floatb_twist_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2IC_gravloadJ_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar2IC_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:37:51
% EndTime: 2020-04-24 20:37:51
% DurationCPUTime: 0.07s
% Computational Cost: add. (33->19), mult. (44->27), div. (3->2), fcn. (19->8), ass. (0->12)
t11 = cos(qJ(2));
t5 = mrSges(3,1) * g(1) + mrSges(3,2) * g(2);
t6 = mrSges(3,1) * g(2) - mrSges(3,2) * g(1);
t9 = sin(qJ(2));
t15 = -t5 * t11 - t6 * t9;
t14 = -qJ(3) + qJ(1);
t13 = t6 * t11 - t5 * t9;
t12 = cos(qJ(1));
t10 = sin(qJ(1));
t8 = m(3) * pkin(2) + mrSges(2,1);
t7 = sin(qJ(2) + t14);
t1 = [(mrSges(2,2) * g(1) - g(2) * t8 + t13) * t12 + t10 * (mrSges(2,2) * g(2) + t8 * g(1) + t15) + ((-pkin(1) * t7 + pkin(2) * sin(t14)) / pkin(1) * (t10 * t15 + t13 * t12) + t9 * ((mrSges(4,1) * g(1) + mrSges(4,2) * g(2)) * sin(qJ(3)) + (-mrSges(4,1) * g(2) + mrSges(4,2) * g(1)) * cos(qJ(3)))) / t7;];
taug = t1(:);
