% Calculate Gravitation load with newton euler on the joints for
% palh2m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:04
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m1IC_gravloadJ_floatb_twist_snew_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1IC_gravloadJ_floatb_twist_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1IC_gravloadJ_floatb_twist_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1IC_gravloadJ_floatb_twist_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1IC_gravloadJ_floatb_twist_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1IC_gravloadJ_floatb_twist_snew_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:04:19
% EndTime: 2020-05-03 01:04:19
% DurationCPUTime: 0.34s
% Computational Cost: add. (527->53), mult. (566->65), div. (0->0), fcn. (398->10), ass. (0->36)
t754 = sin(qJ(2));
t753 = sin(qJ(3));
t758 = cos(qJ(3));
t751 = sin(qJ(5));
t756 = cos(qJ(5));
t742 = pkin(4) * m(6) + t756 * mrSges(6,1) - mrSges(6,2) * t751 + mrSges(5,1);
t749 = m(6) * pkin(6) - mrSges(5,2) + mrSges(6,3);
t752 = sin(qJ(4));
t757 = cos(qJ(4));
t732 = t752 * t742 - t749 * t757;
t794 = mrSges(4,2) + t732;
t733 = t742 * t757 + t749 * t752;
t761 = m(5) + m(6);
t798 = pkin(3) * t761 + mrSges(4,1) + t733;
t727 = t753 * t798 + t758 * t794;
t801 = mrSges(3,2) + t727;
t805 = t754 * t801;
t800 = -t794 * t753 + t798 * t758;
t795 = m(4) + t761;
t725 = t795 * pkin(2) + mrSges(3,1) + t800;
t759 = cos(qJ(2));
t804 = t725 * t759;
t797 = -t732 * t753 + t733 * t758;
t791 = mrSges(2,1) + (m(3) + t795) * pkin(1) + t804;
t755 = sin(qJ(1));
t760 = cos(qJ(1));
t788 = -g(1) * t755 + t760 * g(2);
t787 = g(1) * t760 + t755 * g(2);
t741 = t751 * mrSges(6,1) + mrSges(6,2) * t756 + mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t736 = t754 * g(3) - t787 * t759;
t735 = t759 * g(3) + t787 * t754;
t728 = t732 * t758 + t733 * t753;
t724 = t757 * (t753 * t735 + t758 * t736) + t752 * (t758 * t735 - t753 * t736);
t723 = t728 * t759 + t754 * t797;
t722 = (t754 * t728 - t797 * t759) * g(3);
t1 = [-(-t741 * t760 - t791 * t755) * g(1) - (-t755 * t741 + t791 * t760) * g(2) + t788 * t805; -(-t804 + t805) * g(3) + t722 + t787 * (t725 * t754 + t759 * t801 - t723); -(t727 * t754 - t800 * t759) * g(3) + t722 + t787 * (t727 * t759 + t754 * t800 - t723); mrSges(6,1) * (-t751 * t724 - t756 * t788) - mrSges(6,2) * (t756 * t724 - t751 * t788);];
taug = t1(:);
