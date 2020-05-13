% Implicit kinematic constraints of
% fourbar1turnIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% h [(no of constraints)x1]
%   Implicit constraint equations (e.g. closed loops)
%   In a valid robot configuration qJ this has to be zero

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function h = fourbar1turnIC_kinconstr_impl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_kinconstr_impl_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_kinconstr_impl_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:11
% EndTime: 2020-05-07 11:33:11
% DurationCPUTime: 0.04s
% Computational Cost: add. (11->9), mult. (6->6), div. (0->0), fcn. (6->6), ass. (0->2)
t1 = qJ(2) + qJ(3);
t2 = [cos(qJ(4)) * pkin(4) + pkin(1) + pkin(3) * cos(t1) - cos(qJ(2)) * pkin(2); sin(qJ(4)) * pkin(4) + pkin(3) * sin(t1) - sin(qJ(2)) * pkin(2); pi + qJ(4) - qJ(5) - t1;];
h = t2;
