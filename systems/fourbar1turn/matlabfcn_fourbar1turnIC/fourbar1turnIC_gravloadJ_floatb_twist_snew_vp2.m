% Calculate Gravitation load with newton euler on the joints for
% fourbar1turnIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
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
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnIC_gravloadJ_floatb_twist_snew_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_snew_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:36
% EndTime: 2020-05-07 11:33:36
% DurationCPUTime: 0.11s
% Computational Cost: add. (114->36), mult. (204->52), div. (4->3), fcn. (171->10), ass. (0->23)
t85 = m(4) * pkin(2);
t84 = cos(qJ(1));
t83 = -qJ(4) + qJ(2);
t78 = sin(qJ(1));
t71 = -t84 * g(1) - t78 * g(2);
t77 = sin(qJ(2));
t81 = cos(qJ(2));
t67 = -t81 * g(3) - t77 * t71;
t69 = -t77 * g(3) + t81 * t71;
t76 = sin(qJ(3));
t80 = cos(qJ(3));
t64 = -t80 * t67 + t76 * t69;
t65 = -t76 * t67 - t80 * t69;
t82 = mrSges(4,1) * t64 - mrSges(4,2) * t65;
t79 = cos(qJ(4));
t75 = sin(qJ(4));
t73 = sin(qJ(3) + t83);
t70 = t78 * g(1) - t84 * g(2);
t68 = -t75 * g(3) + t79 * t71;
t66 = -t79 * g(3) - t75 * t71;
t62 = mrSges(4,1) * t70 + mrSges(4,3) * t65;
t61 = -mrSges(4,2) * t70 - mrSges(4,3) * t64;
t1 = [-mrSges(2,2) * t71 + t77 * (-mrSges(3,3) * t67 - t80 * t61 + t76 * t62) + t81 * (mrSges(3,3) * t69 - t76 * t61 - t80 * t62) + (-t75 * t66 + t79 * t68) * mrSges(5,3) + (mrSges(2,1) - t77 * mrSges(3,2) + t81 * (mrSges(3,1) + t85) - t75 * mrSges(5,2) + t79 * mrSges(5,1) + m(5) * pkin(1)) * t70; mrSges(3,1) * t67 - mrSges(3,2) * t69 + (-t73 * t82 + (sin(t83) / pkin(3) * t82 + t76 / pkin(4) * (mrSges(5,1) * t66 - mrSges(5,2) * t68)) * pkin(2)) / t73 + (-t64 * t80 - t65 * t76) * t85 + t82;];
taug = t1(:);
