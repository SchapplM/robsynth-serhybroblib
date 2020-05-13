% Calculate Gravitation load on the joints for
% fourbar2OL
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
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar2OL_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(2,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2OL_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_gravloadJ_floatb_twist_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2OL_gravloadJ_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar2OL_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:19
% EndTime: 2020-04-24 20:32:20
% DurationCPUTime: 0.08s
% Computational Cost: add. (23->14), mult. (40->23), div. (0->0), fcn. (14->6), ass. (0->10)
t4 = mrSges(3,1) * g(1) + mrSges(3,2) * g(2);
t5 = mrSges(3,1) * g(2) - mrSges(3,2) * g(1);
t7 = sin(qJ(2));
t9 = cos(qJ(2));
t12 = -t4 * t9 - t5 * t7;
t11 = -t4 * t7 + t5 * t9;
t10 = cos(qJ(1));
t8 = sin(qJ(1));
t6 = m(3) * pkin(2) + mrSges(2,1);
t1 = [(mrSges(2,2) * g(1) - g(2) * t6 + t11) * t10 + t8 * (mrSges(2,2) * g(2) + t6 * g(1) + t12), t11 * t10 + t8 * t12, (mrSges(4,1) * g(1) + mrSges(4,2) * g(2)) * sin(qJ(3)) + (-mrSges(4,1) * g(2) + mrSges(4,2) * g(1)) * cos(qJ(3)), 0];
taug = t1(:);
